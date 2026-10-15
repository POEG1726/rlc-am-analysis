from __future__ import annotations

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json, math, struct, time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyvisa
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

ROOT = _bench.workspace
SCOPE_ADDR = _bench.scope_address
GEN_ADDR = _bench.generator_address
R_EXT = 5.1  # nominal external resistor; actual DMM measurement not supplied
NPOINTS = 100_000
FS_TARGET = 50e6
COARSE = np.arange(450_000.0, 550_000.1, 2_000.0)

log=[]
def say(s):
    print(s, flush=True); log.append(str(s))

def query_block(inst, command):
    inst.write(command); raw=inst.read_raw(); pos=raw.find(b'#')
    if pos < 0: raise RuntimeError(f'No IEEE binary block for {command!r}: {raw[:60]!r}')
    nd=int(chr(raw[pos+1])); start=pos+2; end=start+nd
    length=int(raw[start:end].decode('ascii')); data=raw[end:end+length]
    if len(data)!=length: raise RuntimeError(f'{command}: expected {length} bytes, got {len(data)}')
    return data

def acquire_channel(scope, ch):
    scope.write(f':WAV:SOUR C{ch}'); scope.write(':WAV:WIDT WORD'); scope.write(':WAV:BYT LSB')
    scope.write(':WAV:INT 1'); scope.write(':WAV:STAR 0'); scope.write(f':WAV:POIN {NPOINTS}')
    pre=query_block(scope, ':WAV:PRE?')
    if len(pre)<188: raise RuntimeError(f'CH{ch}: short preamble {len(pre)}')
    vdiv, voffset, codes_per_div, dt = [struct.unpack_from('<f', pre, p)[0] for p in (156,160,164,176)]
    raw=query_block(scope, ':WAV:DATA?'); codes=np.frombuffer(raw,dtype='<i2')
    if len(codes) < int(0.9 * NPOINTS):
        raise RuntimeError(f'CH{ch}: incomplete waveform block ({len(codes)} of {NPOINTS} points)')
    volts=codes.astype(np.float64)*(vdiv/codes_per_div)-voffset
    return volts, dt, {'vdiv_V_div':vdiv,'voffset_V':voffset,'codes_per_div':codes_per_div,'dt_s':dt,'n':len(volts),'code_min':int(codes.min()),'code_max':int(codes.max())}

def estimate_frequency(x, fs, center):
    # Hann-windowed FFT with log-parabolic interpolation around the local peak.
    y=np.fft.rfft((x-np.mean(x))*np.hanning(len(x)))
    freqs=np.fft.rfftfreq(len(x),1/fs)
    inds=np.flatnonzero((freqs>=center-2500)&(freqs<=center+2500))
    if not len(inds): raise RuntimeError('No bins in frequency search window')
    k=int(inds[np.argmax(np.abs(y[inds]))]); delta=0.0
    if 0<k<len(y)-1:
        a,b,c=np.log(np.maximum(np.abs(y[k-1:k+2]),1e-300))
        den=a-2*b+c
        if den: delta=0.5*(a-c)/den
    return (k+delta)*fs/len(x)

def sine_fit(x, t, f):
    w=2*np.pi*f*t
    A=np.column_stack((np.cos(w),np.sin(w),np.ones(len(t))))
    coeff=np.linalg.lstsq(A,x,rcond=None)[0]
    a,b,dc=coeff
    residual=x-A@coeff
    amp=float(np.hypot(a,b))
    # x=A*cos(wt+phase); positive-frequency complex peak amplitude a - j*b.
    phasor=complex(a,-b)
    return amp, phasor, float(dc), float(np.sqrt(np.mean(residual**2)))

def measure(scope, gen, f_set, label):
    gen.write(f'C2:BSWV FRQ,{f_set:.6f}HZ')
    time.sleep(0.10)
    last_error = None
    for attempt in range(3):
        # The SDS824X HD occasionally returns an empty IEEE block if stopped
        # before its auto-triggered acquisition has completed. Allow time for
        # the 100 kpoint record to become available before stopping it.
        scope.write(':TRIG:RUN'); time.sleep(0.50); scope.write(':TRIG:STOP'); time.sleep(0.01)
        try:
            x1,dt1,i1=acquire_channel(scope,1); x2,dt2,i2=acquire_channel(scope,2)
            last_error = None
            break
        except RuntimeError as e:
            last_error = e
            if attempt == 2: raise
            time.sleep(0.10)
    if last_error is not None: raise last_error
    if not np.isclose(dt1,dt2,rtol=0,atol=1e-15): raise RuntimeError(f'CH sample interval mismatch {dt1} vs {dt2}')
    n=min(len(x1),len(x2)); x1=x1[:n]; x2=x2[:n]; fs=1/dt1; t=np.arange(n)/fs
    fest=estimate_frequency(x1,fs,f_set)
    ain,cin,dcin,rmsin=sine_fit(x1,t,fest); aout,cout,dcout,rmsout=sine_fit(x2,t,fest)
    H=cout/cin; Z=R_EXT/H
    row={'label':label,'frequency_set_Hz':f_set,'frequency_measured_Hz':fest,'sample_rate_Hz':fs,'dt_s':dt1,'n_points':n,'record_duration_s':n*dt1,
         'vin_Vpk':ain,'vout_Vpk':aout,'H_real':H.real,'H_imag':H.imag,'H_mag':abs(H),'H_phase_deg':np.angle(H,deg=True),
         'Z_real_ohm':Z.real,'Z_imag_ohm':Z.imag,'Z_mag_ohm':abs(Z),'THD_residual_CH1_Vrms':rmsin,'THD_residual_CH2_Vrms':rmsout,
         'vin_dc_V':dcin,'vout_dc_V':dcout,'CH1_vdiv_V_div':i1['vdiv_V_div'],'CH2_vdiv_V_div':i2['vdiv_V_div'],
         'CH1_code_min':i1['code_min'],'CH1_code_max':i1['code_max'],'CH2_code_min':i2['code_min'],'CH2_code_max':i2['code_max']}
    return row

def fit_complex(df):
    # R_s, L, C, and a common channel timing offset. Fit a symmetric region around measured peak.
    d=df[np.isfinite(df.H_real)&np.isfinite(df.H_imag)].copy()
    fpk=float(d.loc[d.H_mag.idxmax(),'frequency_measured_Hz'])
    d=d[np.abs(d.frequency_measured_Hz-fpk)<=30_000]
    f=d.frequency_measured_Hz.to_numpy(); hm=d.H_real.to_numpy()+1j*d.H_imag.to_numpy()
    def residual(p):
        rs,L,C,tau=p
        z=rs+1j*(2*np.pi*f*L-1/(2*np.pi*f*C))
        h=R_EXT/z*np.exp(-1j*2*np.pi*f*tau)
        return np.r_[((h-hm).real),((h-hm).imag)]
    p0=[7.7,22e-6,4.7e-9,0.0]
    fit=least_squares(residual,p0,bounds=([R_EXT*0.7,1e-6,.5e-9,-2e-6],[R_EXT*10,100e-6,20e-9,2e-6]),x_scale=[5,20e-6,5e-9,1e-9],max_nfev=10000)
    rs,L,C,tau=fit.x; residuals=residual(fit.x)
    pred=R_EXT/(rs+1j*(2*np.pi*f*L-1/(2*np.pi*f*C)))*np.exp(-1j*2*np.pi*f*tau)
    rmse=float(np.sqrt(np.mean(np.abs(pred-hm)**2)))
    return {'Rs_fit_ohm':rs,'L_fit_H':L,'C_fit_F':C,'channel_delay_fit_s':tau,'complex_H_RMSE':rmse,'fit_n':len(f),'fit_success':bool(fit.success),'fit_message':fit.message,'fit_frequency_min_Hz':float(f.min()),'fit_frequency_max_Hz':float(f.max())}

def crossing(f,h,idx,step):
    target=h[idx]/np.sqrt(2)
    lo=idx
    while lo>0 and h[lo]>=target: lo-=1
    hi=idx
    while hi<len(h)-1 and h[hi]>=target: hi+=1
    if lo==idx or hi==idx: return math.nan,math.nan
    fl=np.interp(target,[h[lo],h[lo+1]],[f[lo],f[lo+1]])
    fu=np.interp(target,[h[hi],h[hi-1]],[f[hi],f[hi-1]])
    return float(fl),float(fu)

def main():
    rm=pyvisa.ResourceManager(); sc=rm.open_resource(SCOPE_ADDR); gn=rm.open_resource(GEN_ADDR)
    sc.timeout=120000; gn.timeout=8000; sc.chunk_size=4*1024*1024
    rows=[]; initial={}
    try:
        say(f'UTC start: {datetime.now(timezone.utc).isoformat()}')
        say('Scope: '+sc.query('*IDN?').strip()); say('Generator: '+gn.query('*IDN?').strip())
        initial['generator_bswv']=gn.query('C2:BSWV?').strip(); initial['generator_mdwv_before']=gn.query('C2:MDWV?').strip(); initial['generator_outp']=gn.query('C2:OUTP?').strip()
        say('CH2 before: '+initial['generator_bswv']); say('AM before: '+initial['generator_mdwv_before']); say('CH2 output: '+initial['generator_outp'])
        # User explicitly requested AM OFF and pure sine; set CH2 only.
        gn.write('C2:MDWV STATE,OFF'); time.sleep(.15)
        if 'STATE,OFF' not in gn.query('C2:MDWV?').upper(): raise RuntimeError('AM did not turn OFF on generator CH2')
        gn.write('C2:BSWV WVTP,SINE,AMP,5V,OFST,0V')
        gn.write('C2:OUTP ON,LOAD,HZ')
        say('AM after: '+gn.query('C2:MDWV?').strip()); say('CH2 sine config: '+gn.query('C2:BSWV?').strip())
        # Short high-rate diagnostic record; no transfer decimation.
        sc.write(':TIM:SCAL 2E-4'); time.sleep(.2)
        say(f'Scope settings: points={sc.query(":ACQ:POIN?").strip()}, sample_rate={sc.query(":ACQ:SRAT?").strip()}, timebase={sc.query(":TIM:SCAL?").strip()}')
        for ch in (1,2):
            say(f'CH{ch} config: scale={sc.query(f":CHAN{ch}:SCAL?").strip()} V/div, impedance={sc.query(f":CHAN{ch}:IMP?").strip()}, probe={sc.query(f":CHAN{ch}:PROB?").strip()}, coupling={sc.query(f":CHAN{ch}:COUP?").strip()}, bandwidth={sc.query(f":CHAN{ch}:BWL?").strip()}')
        for f in COARSE:
            row=measure(sc,gn,float(f),'coarse'); rows.append(row)
            say(f'coarse {f/1000:.1f} kHz: measured {row["frequency_measured_Hz"]/1000:.4f} kHz H={row["H_mag"]:.6f} phase={row["H_phase_deg"]:+.3f} deg Vin={row["vin_Vpk"]:.4g}Vpk Vout={row["vout_Vpk"]:.4g}Vpk')
        coarse_df=pd.DataFrame(rows); ix=int(coarse_df.H_mag.idxmax()); fpk=float(coarse_df.loc[ix,'frequency_measured_Hz']); hp=float(coarse_df.loc[ix,'H_mag'])
        # Estimate rough bandwidth from coarse -3dB crossings, then oversample with 500 Hz points.
        fc=coarse_df.frequency_measured_Hz.to_numpy(); hc=coarse_df.H_mag.to_numpy(); flo,fhi=crossing(fc,hc,ix,2000)
        rough_bw=fhi-flo if np.isfinite(flo+fhi) else 40_000
        span=min(45_000.0,max(12_000.0,1.2*rough_bw))
        fine_lo=max(450_000.0,fpk-span); fine_hi=min(550_000.0,fpk+span)
        fine=np.arange(math.ceil(fine_lo/500)*500, fine_hi+250,500,dtype=float)
        say(f'Coarse peak at {fpk:.2f} Hz; coarse -3 dB estimates {flo:.2f}, {fhi:.2f} Hz; fine sweep {fine_lo:.1f}..{fine_hi:.1f} Hz step 500 Hz ({len(fine)} points)')
        for f in fine:
            row=measure(sc,gn,float(f),'fine'); rows.append(row)
            if len(rows)%20==0: say(f'fine progress {sum(r["label"]=="fine" for r in rows)}/{len(fine)}; last {f/1000:.2f} kHz H={row["H_mag"]:.6f}')
        # Three repeat captures at approximate resonance to quantify short-term variation.
        for rep in range(3):
            row=measure(sc,gn,fpk,f'repeat_{rep+1}'); rows.append(row)
            say(f'repeat {rep+1}: f={row["frequency_measured_Hz"]:.2f} Hz H={row["H_mag"]:.7f} phase={row["H_phase_deg"]:+.4f} deg')
        df=pd.DataFrame(rows); df.to_csv(ROOT/'resonance_scan.csv',index=False)
        # Main response curve is fine scan; direct interpolate -3dB points.
        fine_df=df[df.label=='fine'].sort_values('frequency_measured_Hz'); f=fine_df.frequency_measured_Hz.to_numpy(); h=fine_df.H_mag.to_numpy(); ip=int(np.argmax(h)); f0_direct=float(f[ip]); hmax=float(h[ip]); fl,fu=crossing(f,h,ip,500); bw=fu-fl; q=f0_direct/bw
        fit=fit_complex(df[df.label.isin(['coarse','fine'])]); rs=fit['Rs_fit_ohm']; L=fit['L_fit_H']; C=fit['C_fit_F']
        f0_model=1/(2*np.pi*np.sqrt(L*C)); bw_model=rs/(2*np.pi*L); q_model=f0_model/bw_model
        # Direct Hmax-derived resistance from known Rext, with fitted AC loss excluding Rext.
        rs_direct=R_EXT/hmax
        # Diagnostics and model columns
        zpred=rs+1j*(2*np.pi*fine_df.frequency_measured_Hz.to_numpy()*L-1/(2*np.pi*fine_df.frequency_measured_Hz.to_numpy()*C))
        tau=fit['channel_delay_fit_s']; hpred=R_EXT/zpred*np.exp(-1j*2*np.pi*fine_df.frequency_measured_Hz.to_numpy()*tau)
        fine_df=fine_df.copy(); fine_df['H_model_mag']=np.abs(hpred); fine_df['H_model_phase_deg']=np.angle(hpred,deg=True)
        # Plot magnitude response and -3 dB levels.
        fig,ax=plt.subplots(figsize=(10,6)); ax.plot(f/1000,h,'o',ms=3,label='Measured |Vout/Vin|'); ax.plot(f/1000,np.abs(hpred),'-',lw=1.6,label='Complex RLC fit')
        ax.axvline(f0_direct/1000,color='k',ls='--',alpha=.55,label=f'Measured peak {f0_direct/1000:.3f} kHz'); ax.axhline(hmax/np.sqrt(2),color='gray',ls=':')
        if np.isfinite(fl): ax.axvline(fl/1000,color='gray',ls=':'); ax.axvline(fu/1000,color='gray',ls=':',label=f'-3 dB: {fl/1000:.3f}, {fu/1000:.3f} kHz')
        ax.set(xlabel='Measured frequency (kHz)',ylabel='|Vout/Vin|',title='Series RLC resistor transfer response'); ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'resonance_curve.png',dpi=180); plt.close(fig)
        fig,ax=plt.subplots(figsize=(10,6)); ax.plot(df.frequency_measured_Hz/1000,df.Z_real_ohm,'o',ms=3,label='Measured Re(Z)'); ax.plot(df.frequency_measured_Hz/1000,df.Z_imag_ohm,'o',ms=3,label='Measured Im(Z)')
        fplot=np.linspace(450e3,550e3,500); zmod=rs+1j*(2*np.pi*fplot*L-1/(2*np.pi*fplot*C)); ax.plot(fplot/1000,zmod.real,'-',lw=1.3,label='Fit Re(Z)'); ax.plot(fplot/1000,zmod.imag,'-',lw=1.3,label='Fit Im(Z)')
        ax.axhline(0,color='k',lw=.6); ax.set(xlabel='Measured frequency (kHz)',ylabel='Impedance (ohm)',title='Complex impedance derived from two-channel phasors'); ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'impedance_real_imag.png',dpi=180); plt.close(fig)
        # Repeated point statistics.
        reps=df[df.label.str.startswith('repeat_')]
        repeat_h_sd=float(reps.H_mag.std(ddof=1)); repeat_h_mean=float(reps.H_mag.mean()); repeat_phase_sd=float(reps.H_phase_deg.std(ddof=1))
        params={'f0_measured_peak_Hz':f0_direct,'Hmax_measured':hmax,'f_lower_3dB_Hz':fl,'f_upper_3dB_Hz':fu,'BW_measured_Hz':bw,'Q_measured':q,
                'Rs_eff_from_Hmax_ohm':rs_direct,'R_external_nominal_ohm':R_EXT,'series_loss_excluding_external_fit_ohm':rs-R_EXT,
                **fit,'f0_from_fit_Hz':f0_model,'BW_from_fit_Hz':bw_model,'Q_from_fit':q_model,
                'L_eff_fit_H':L,'C_eff_fit_F':C,'L_eff_fit_uH':L*1e6,'C_eff_fit_nF':C*1e9,
                'Rs_total_fit_ohm':rs,'repeat_H_mean':repeat_h_mean,'repeat_H_sd':repeat_h_sd,'repeat_H_sd_percent':100*repeat_h_sd/repeat_h_mean,'repeat_phase_sd_deg':repeat_phase_sd,
                'nominal_L_H':22e-6,'nominal_C_F':4.7e-9,'nominal_f0_Hz':1/(2*np.pi*np.sqrt(22e-6*4.7e-9)),
                'scan_utc':datetime.now(timezone.utc).isoformat(),'scope_id':sc.query('*IDN?').strip(),'generator_id':gn.query('*IDN?').strip()}
        (ROOT/'characterization_parameters.json').write_text(json.dumps(params,indent=2),encoding='utf-8')
        say('PARAMETERS '+json.dumps(params,indent=2))
        # Re-evaluate old 500k +/- 1k using measured response curve; compare to user reported result.
        def interp_h(freq): return float(np.interp(freq,f,h))
        predicted_tlo=interp_h(500_000-1000)/interp_h(500_000); predicted_thi=interp_h(500_000+1000)/interp_h(500_000)
        say(f'500 kHz detuning check: measured/interpolated T(-1 kHz)={predicted_tlo:.6f}, T(+1 kHz)={predicted_thi:.6f}, mean={(predicted_tlo+predicted_thi)/2:.6f}')
        gn.write(f'C2:BSWV FRQ,{f0_direct:.6f}HZ'); time.sleep(.15)
        initial['mdwv_final']=gn.query('C2:MDWV?').strip(); initial['bswv_final']=gn.query('C2:BSWV?').strip(); initial['output_final']=gn.query('C2:OUTP?').strip()
        say('Final generator state: '+initial['bswv_final']+'; '+initial['mdwv_final']+'; '+initial['output_final'])
        say('UTC end: '+datetime.now(timezone.utc).isoformat())
        (ROOT/'characterization_log.txt').write_text('\n'.join(log)+'\n',encoding='utf-8')
    finally:
        # Restore requested long-record scope acquisition geometry and leave CH2 at measured resonance.
        try:
            sc.write(':TRIG:RUN'); sc.write(':TIM:SCAL 5E-2')
            gn.write('C2:MDWV STATE,OFF')
            sc.close(); gn.close(); rm.close()
        except Exception: pass

if __name__=='__main__': main()
