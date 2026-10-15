"""One-point 1 kHz / 50% AM diagnostic at the characterized low-drive RLC level.
Based on diagnose_am_1khz.py; preserves a separate output set and restores the
previous 5 Vpp / 493.9 kHz / AM-OFF bench configuration on exit.
"""

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json, struct, time
from pathlib import Path
import numpy as np, pandas as pd, pyvisa
import matplotlib.pyplot as plt
ROOT=_bench.workspace
SCOPE=_bench.scope_address
GEN=_bench.generator_address
N=5_000_000; FM=1000.0; DEPTH=50.0; AMP=0.20; FC=498337.4
# Empirical relative range factors to the 200 mV/div reference from the matched-sine calibration.
cal=json.loads((ROOT/'vertical_scale_corrections.json').read_text(encoding='utf-8'))['correction_factors_to_reference']
def corr(ch,scale):
    d=cal[f'CH{ch}']; return float(d[str(float(scale))] if str(float(scale)) in d else d[f'{float(scale):g}'])
log=[]
def say(s): print(s,flush=True);log.append(str(s))
def block(inst,cmd):
    inst.write(cmd); raw=inst.read_raw(); p=raw.find(b'#')
    if p<0: raise RuntimeError(f'No binary block for {cmd}: {raw[:60]!r}')
    nd=int(chr(raw[p+1])); nbytes=int(raw[p+2:p+2+nd]); data=raw[p+2+nd:p+2+nd+nbytes]
    if len(data)!=nbytes: raise RuntimeError(f'Truncated {cmd}: {len(data)}/{nbytes}')
    return data
def capture(sc,ch):
    sc.write(f':WAV:SOUR C{ch}');sc.write(':WAV:WIDT WORD');sc.write(':WAV:BYT LSB');sc.write(':WAV:INT 1');sc.write(':WAV:STAR 0');sc.write(f':WAV:POIN {N}')
    pre=block(sc,':WAV:PRE?');vdiv,offset,cpd,dt=[struct.unpack_from('<f',pre,p)[0] for p in (156,160,164,176)]
    codes=np.frombuffer(block(sc,':WAV:DATA?'),dtype='<i2');v=codes.astype(np.float64)*(vdiv/cpd)-offset
    return v,dt,dict(vdiv_V_div=vdiv,offset_V=offset,codes_per_div=cpd,n=len(codes),min_code=int(codes.min()),max_code=int(codes.max()))
def peak(freq,amp,target,width=100):
    ids=np.flatnonzero(np.abs(freq-target)<=width)
    if not len(ids): raise RuntimeError(f'No FFT bins near {target:g} Hz')
    k=int(ids[np.argmax(amp[ids])]);delta=0.0
    if 0<k<len(amp)-1:
        a,b,c=np.log(np.maximum(amp[k-1:k+2],1e-300));den=a-2*b+c
        if den: delta=.5*(a-c)/den
    return float((k+delta)*(freq[1]-freq[0])),float(amp[k]),k
rm=sc=gn=None
try:
    rm=pyvisa.ResourceManager();sc=rm.open_resource(SCOPE);gn=rm.open_resource(GEN)
    sc.timeout=120000;gn.timeout=10000;sc.chunk_size=16*1024*1024
    say('Scope '+sc.query('*IDN?').strip());say('Generator '+gn.query('*IDN?').strip())
    say('PRE generator '+gn.query('C2:BSWV?').strip()+'; '+gn.query('C2:MDWV?').strip()+'; '+gn.query('C2:OUTP?').strip())
    say(f'PRE scope timebase={sc.query(":TIM:SCAL?").strip()}, CH1={sc.query(":CHAN1:SCAL?").strip()} V/div, CH2={sc.query(":CHAN2:SCAL?").strip()} V/div, points={sc.query(":ACQ:POIN?").strip()}, rate={sc.query(":ACQ:SRAT?").strip()}')
    # At the AM envelope peak, measured Vin is expected to be under 4 divisions at 10 mV/div;
    # CH2 is kept at 5 mV/div to avoid the characterized 2 mV/div gain error.
    sc.write(':TIM:SCAL 5E-2');sc.write(':CHAN1:SCAL 2E-2');sc.write(':CHAN2:SCAL 1E-2')
    gn.write('C2:MDWV STATE,OFF');time.sleep(.15)
    gn.write(f'C2:BSWV WVTP,SINE,FRQ,{FC:.1f}HZ,AMP,{AMP:.4g}V,OFST,0V');time.sleep(.15)
    gn.write(f'C2:MDWV STATE,ON,AM,MDSP,SINE,SRC,INT,FRQ,{FM:.1f}HZ,DEPTH,{DEPTH:g}');time.sleep(.3)
    md=gn.query('C2:MDWV?').strip(); bswv=gn.query('C2:BSWV?').strip(); outp=gn.query('C2:OUTP?').strip()
    say('AM '+md);say('Carrier '+bswv);say('Output '+outp)
    if 'STATE,ON' not in md.upper() or 'FRQ,1000HZ' not in md.upper() or 'DEPTH,50' not in md.upper(): raise RuntimeError('AM setting readback does not match 1 kHz / 50% request')
    say(f'Capture config points={sc.query(":ACQ:POIN?").strip()} rate={sc.query(":ACQ:SRAT?").strip()} timebase={sc.query(":TIM:SCAL?").strip()} CH1/CH2 scales={sc.query(":CHAN1:SCAL?").strip()}/{sc.query(":CHAN2:SCAL?").strip()} V/div')
    sc.write(':TRIG:RUN');time.sleep(.6);sc.write(':TRIG:STOP');time.sleep(.05)
    vin,dt1,m1=capture(sc,1);vout,dt2,m2=capture(sc,2)
    if len(vin)!=N or len(vout)!=N or not np.isclose(dt1,dt2): raise RuntimeError(f'Bad paired waveform: {len(vin)}, {len(vout)}, dt={dt1}/{dt2}')
    if m1['min_code']<=-32768 or m1['max_code']>=32767 or m2['min_code']<=-32768 or m2['max_code']>=32767: raise RuntimeError('ADC clipping detected; refusing to analyze saturated capture')
    fs=1/dt1;n=len(vin);df=fs/n;duration=n/fs
    say(f'Acquired paired records: N={n}, dt={dt1:.12g}s, fs={fs:.6f}Sa/s, duration={duration:.6f}s, FFT bin={df:.6f}Hz')
    # Existing diagnostic method: Hann window, independent local peak interpolation on each channel.
    w=np.hanning(n);cg=w.mean();freq=np.fft.rfftfreq(n,dt1)
    def spectrum(x):
        x=x-np.mean(x);a=2*np.abs(np.fft.rfft(x*w))/(n*cg);a[0]/=2;return a
    ain=spectrum(vin);aout=spectrum(vout)
    rows=[]
    for label,target in [('LSB',FC-FM),('carrier',FC),('USB',FC+FM)]:
        fi,vi,ki=peak(freq,ain,target);fo,vo,ko=peak(freq,aout,target)
        f1=corr(1,m1['vdiv_V_div']);f2=corr(2,m2['vdiv_V_div'])
        rows.append({'component':label,'target_Hz':target,'Vin_fft_Vpk_raw':vi,'Vin_fft_Vpk_corrected':vi*f1,'Vin_fft_peak_Hz':fi,'Vout_fft_Vpk_raw':vo,'Vout_fft_Vpk_corrected':vo*f2,'Vout_fft_peak_Hz':fo,'H_raw':vo/vi,'H_range_corrected':vo*f2/(vi*f1),'in_bin':ki,'out_bin':ko,'same_bin_H_raw':aout[int(round(target/df))]/ain[int(round(target/df))]})
        say(f'{label}: target={target:.1f}Hz Vin={vi*f1:.7g}Vpk(corrected) at {fi:.2f}Hz; Vout={vo*f2:.7g}Vpk(corrected) at {fo:.2f}Hz; H={vo*f2/(vi*f1):.7f}; bins={ki}/{ko}')
    h=np.array([r['H_range_corrected'] for r in rows]);tl,tc,tu=h
    tlo,thi=tl/tc,tu/tc;tmean=(tlo+thi)/2
    say(f'Normalized transfer: T_lower={tlo:.7f}, T_upper={thi:.7f}, T_mean={tmean:.7f}')
    say(f'Raw scope range CH1={vin.min():.6g}..{vin.max():.6g}V, CH2={vout.min():.6g}..{vout.max():.6g}V; dBc results are unaffected by the relative range gain factors')
    say('CH1 acquisition metadata '+json.dumps(m1));say('CH2 acquisition metadata '+json.dumps(m2))
    dfout=pd.DataFrame(rows);dfout['amplitude_programmed_Vpp']=AMP;dfout['carrier_programmed_Hz']=FC;dfout['modulation_frequency_Hz']=FM;dfout['modulation_depth_percent']=DEPTH;dfout['record_points']=n;dfout['sample_rate_Hz']=fs;dfout['fft_bin_Hz']=df;dfout.to_csv(ROOT/'am_1khz_low_drive.csv',index=False)
    np.savez_compressed(ROOT/'am_1khz_low_drive.npz',vin=vin.astype(np.float32),vout=vout.astype(np.float32),dt_s=dt1,fc_set_Hz=FC,fm_set_Hz=FM,amplitude_set_Vpp=AMP,depth_percent=DEPTH)
    pd.DataFrame([{'carrier_programmed_Hz':FC,'modulation_frequency_Hz':FM,'depth_percent':DEPTH,'amplitude_programmed_Vpp':AMP,'T_lower':tlo,'T_upper':thi,'T_mean':tmean,'CH1_scale_V_div':m1['vdiv_V_div'],'CH2_scale_V_div':m2['vdiv_V_div'],'CH1_range_factor':f1,'CH2_range_factor':f2,'points':n,'sample_rate_Hz':fs,'duration_s':duration,'FFT_bin_Hz':df,'Vin_raw_min_V':float(vin.min()),'Vin_raw_max_V':float(vin.max()),'Vout_raw_min_V':float(vout.min()),'Vout_raw_max_V':float(vout.max())}]).to_csv(ROOT/'am_1khz_low_drive_summary.csv',index=False)
    # Focused spectral plot around the triplet.
    mask=(freq>=FC-5000)&(freq<=FC+5000)
    fig,ax=plt.subplots(figsize=(10,5.5));ax.plot((freq[mask]-FC)/1000,20*np.log10(np.maximum(ain[mask],1e-12)),label='CH1 Vin');ax.plot((freq[mask]-FC)/1000,20*np.log10(np.maximum(aout[mask],1e-12)),label='CH2 Vout')
    for x,label in [(-1,'LSB'),(0,'carrier'),(1,'USB')]:ax.axvline(x,ls='--',lw=.8,alpha=.5,label=label)
    ax.set(xlabel='Frequency relative to carrier (kHz)',ylabel='Amplitude (dBV peak, raw scope scaling)',title='Low-drive 1 kHz AM spectrum: 0.2 Vpp, 50% depth');ax.grid(alpha=.25);ax.legend(ncol=2);fig.tight_layout();fig.savefig(ROOT/'am_1khz_low_drive_spectrum.png',dpi=170);plt.close(fig)
    (ROOT/'am_1khz_low_drive.log').write_text('\n'.join(log)+'\n',encoding='utf-8')
finally:
    if sc is not None:
        try: sc.write(':TRIG:RUN');sc.write(':TIM:SCAL 5E-2');sc.write(':CHAN1:SCAL 2E-1');sc.write(':CHAN2:SCAL 1E-1');sc.close()
        except Exception as e: say(f'Scope restore warning: {e!r}')
    if gn is not None:
        try:
            gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,FRQ,493900HZ,AMP,5V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ')
            say('FINAL generator '+gn.query('C2:BSWV?').strip()+'; '+gn.query('C2:MDWV?').strip()+'; '+gn.query('C2:OUTP?').strip());gn.close()
        except Exception as e: say(f'Generator restore warning: {e!r}')
    if rm is not None: rm.close()
    if log: (ROOT/'am_1khz_low_drive.log').write_text('\n'.join(log)+'\n',encoding='utf-8')
