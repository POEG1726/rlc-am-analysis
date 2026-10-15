"""Low-drive AM sideband SNR and three-capture repeatability sweep.
Reuses the validated 5 Mpoint native waveform/FFT method from the existing AM diagnostic.
"""

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json, math, re, struct, time
from pathlib import Path
import numpy as np, pandas as pd, pyvisa
import matplotlib.pyplot as plt
ROOT=_bench.workspace
SCOPE=_bench.scope_address;GEN=_bench.generator_address
FM_LIST=[1000,2000,3000,5000,8000,10000,12000,15000,18000,20000,22000,25000,30000,35000,40000]
FC=498337.4;AMP=0.2;DEPTH=50.0;N=5_000_000;REPS=3;R_EXT=5.1
SCALES={1:0.02,2:0.01};TIMEBASE=0.05
REPEAT_CSV=ROOT/'low_drive_snr_repeats.csv';SUMMARY_CSV=ROOT/'low_drive_snr_sweep.csv';LOG_PATH=ROOT/'low_drive_snr_sweep.log'
cal=json.loads((ROOT/'vertical_scale_corrections.json').read_text(encoding='utf-8'))['correction_factors_to_reference']
def scale_factor(ch,scale):
    return next(float(v) for k,v in cal[f'CH{ch}'].items() if np.isclose(float(k),float(scale),rtol=0,atol=1e-9))
FACTORS={ch:scale_factor(ch,SCALES[ch]) for ch in (1,2)}
log=[]
def say(s): print(s,flush=True);log.append(str(s));LOG_PATH.write_text('\n'.join(log)+'\n',encoding='utf-8')
def query_block(inst,command):
    inst.write(command);raw=inst.read_raw();p=raw.find(b'#')
    if p<0:raise RuntimeError(f'No IEEE block for {command}: {raw[:60]!r}')
    nd=int(chr(raw[p+1]));start=p+2;end=start+nd;size=int(raw[start:end].decode('ascii'));data=raw[end:end+size]
    if len(data)!=size:raise RuntimeError(f'{command}: got {len(data)} of {size} bytes')
    return data
def capture(scope,ch):
    scope.write(f':WAV:SOUR C{ch}');scope.write(':WAV:WIDT WORD');scope.write(':WAV:BYT LSB');scope.write(':WAV:INT 1');scope.write(':WAV:STAR 0');scope.write(f':WAV:POIN {N}')
    pre=query_block(scope,':WAV:PRE?')
    vdiv,offset,cpd,dt=[struct.unpack_from('<f',pre,p)[0] for p in (156,160,164,176)]
    codes=np.frombuffer(query_block(scope,':WAV:DATA?'),dtype='<i2')
    volts=codes.astype(np.float64)*(vdiv/cpd)-offset
    return volts,dt,{'vdiv_V_div':vdiv,'offset_V':offset,'codes_per_div':cpd,'points':len(codes),'min_code':int(codes.min()),'max_code':int(codes.max())}
def spectrum(x,dt):
    n=len(x);window=np.hanning(n);cg=window.mean();z=np.fft.rfft((x-np.mean(x))*window)
    a=2*np.abs(z)/(n*cg);a[0]/=2;freq=np.fft.rfftfreq(n,dt)
    return freq,a
def line_and_noise(freq,amp,target):
    # Local line peak within +/-100 Hz; local floor uses bins 100..450 Hz on either side.
    ids=np.flatnonzero(np.abs(freq-target)<=100.0)
    if not len(ids):raise RuntimeError(f'No FFT samples near {target} Hz')
    k=int(ids[np.argmax(amp[ids])]);line=float(amp[k]);fline=float(freq[k])
    mask=(np.abs(freq-target)>=100.0)&(np.abs(freq-target)<=450.0)
    vals=amp[mask].copy();fnoise=freq[mask]
    if len(vals)<20:raise RuntimeError('Too few local noise bins')
    # Remove isolated narrow spurs robustly: reject a candidate >12 dB above the local median
    # and guard +/-20 Hz around each rejected bin, then take the median remaining bin amplitude.
    med=float(np.median(vals));spur=vals>med*(10**(12/20))
    keep=np.ones(len(vals),dtype=bool)
    for fspur in fnoise[spur]:keep &= np.abs(fnoise-fspur)>20.0
    floor=float(np.median(vals[keep])) if keep.sum()>=10 else med
    snr=20*math.log10(max(line,1e-30)/max(floor,1e-30))
    return fline,line,floor,snr,int(k),int(keep.sum()),int((~keep).sum())
def analyze_channel(freq,a,ch,fm):
    targets={'LSB':FC-fm,'carrier':FC,'USB':FC+fm};out={}
    for name,target in targets.items():
        fl,amp,linefloor,snr,k,nnoise,nspur=line_and_noise(freq,a,target)
        factor=FACTORS[ch]
        out[name]={'frequency':fl,'amp_raw':amp,'amp_corrected':amp*factor,'noise_raw':linefloor,'noise_corrected':linefloor*factor,'snr':snr,'bin':k,'noise_bins':nnoise,'excluded_spur_bins':nspur}
    return out
# Resume safely from per-repetition checkpoints if this script is interrupted.
if REPEAT_CSV.exists():
    prior=pd.read_csv(REPEAT_CSV)
    if len(prior) and not ((prior.amplitude_programmed_Vpp==AMP).all() and (prior.carrier_programmed_Hz==FC).all() and (prior.depth_percent==DEPTH).all() and (prior.record_points==N).all()):
        raise RuntimeError('Existing checkpoint has different conditions; preserve it before starting a new run.')
    done={(int(r.fm_Hz),int(r.repeat)) for r in prior.itertuples()} if len(prior) else set()
else:
    prior=pd.DataFrame();done=set()
rm=sc=gn=None
try:
    rm=pyvisa.ResourceManager();sc=rm.open_resource(SCOPE);gn=rm.open_resource(GEN);sc.timeout=120000;gn.timeout=10000;sc.chunk_size=16*1024*1024
    say('Scope '+sc.query('*IDN?').strip());say('Generator '+gn.query('*IDN?').strip())
    say('Initial generator: '+gn.query('C2:BSWV?').strip()+'; '+gn.query('C2:MDWV?').strip()+'; '+gn.query('C2:OUTP?').strip())
    say(f'Initial scope: timebase={sc.query(":TIM:SCAL?").strip()}, sample_rate={sc.query(":ACQ:SRAT?").strip()}, points={sc.query(":ACQ:POIN?").strip()}, CH1={sc.query(":CHAN1:SCAL?").strip()}, CH2={sc.query(":CHAN2:SCAL?").strip()}')
    gn.write('C2:MDWV STATE,OFF');gn.write(f'C2:BSWV WVTP,SINE,FRQ,{FC:.1f}HZ,AMP,{AMP:.4g}V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ')
    sc.write(f':TIM:SCAL {TIMEBASE:.9g}');sc.write(f':CHAN1:SCAL {SCALES[1]:.9g}');sc.write(f':CHAN2:SCAL {SCALES[2]:.9g}')
    time.sleep(.2)
    say(f'Fixed capture conditions: fc={FC:.1f}Hz, amp={AMP:g}Vpp, AM depth={DEPTH:g}%, scope={sc.query(":ACQ:SRAT?").strip()}Sa/s, {sc.query(":ACQ:POIN?").strip()}points, scales CH1/CH2={sc.query(":CHAN1:SCAL?").strip()}/{sc.query(":CHAN2:SCAL?").strip()}V/div')
    for fm in FM_LIST:
        gn.write(f'C2:MDWV STATE,ON,AM,MDSP,SINE,SRC,INT,FRQ,{fm}HZ,DEPTH,{DEPTH:g}');time.sleep(.25)
        md=gn.query('C2:MDWV?').strip()
        if 'STATE,ON' not in md.upper() or f'FRQ,{fm}HZ' not in md.upper() or 'DEPTH,50' not in md.upper():raise RuntimeError(f'AM readback mismatch at {fm}Hz: {md}')
        for rep in range(1,REPS+1):
            if (fm,rep) in done:continue
            sc.write(':TRIG:RUN');time.sleep(.65);sc.write(':TRIG:STOP');time.sleep(.05)
            vin,dt1,m1=capture(sc,1);vout,dt2,m2=capture(sc,2)
            if len(vin)!=N or len(vout)!=N or not np.isclose(dt1,dt2,rtol=0,atol=1e-15):raise RuntimeError(f'Bad paired record at {fm} Hz rep {rep}')
            clip=(m1['min_code']<=-32768 or m1['max_code']>=32767 or m2['min_code']<=-32768 or m2['max_code']>=32767)
            if clip:raise RuntimeError(f'ADC clipping at {fm} Hz rep {rep}; refusing analysis. codes={m1}/{m2}')
            freq,ain=spectrum(vin,dt1);_,aout=spectrum(vout,dt2)
            ch1=analyze_channel(freq,ain,1,fm);ch2=analyze_channel(freq,aout,2,fm)
            row={'fm_Hz':fm,'repeat':rep,'amplitude_programmed_Vpp':AMP,'carrier_programmed_Hz':FC,'depth_percent':DEPTH,'record_points':N,'sample_rate_Hz':1/dt1,'record_duration_s':N*dt1,'fft_bin_Hz':1/(N*dt1),'CH1_scale_V_div':m1['vdiv_V_div'],'CH2_scale_V_div':m2['vdiv_V_div'],'CH1_scale_correction':FACTORS[1],'CH2_scale_correction':FACTORS[2],'CH1_clipped':False,'CH2_clipped':False,'CH1_code_min':m1['min_code'],'CH1_code_max':m1['max_code'],'CH2_code_min':m2['min_code'],'CH2_code_max':m2['max_code']}
            for ch,data in ((1,ch1),(2,ch2)):
                for name,x in data.items():
                    prefix=f'CH{ch}_{name}'
                    for key,value in x.items():row[f'{prefix}_{key}']=value
            for component in ('LSB','carrier','USB'):
                row[f'{component}_same_bin_CH1_CH2']=ch1[component]['bin']==ch2[component]['bin']
            hl=ch2['LSB']['amp_corrected']/ch1['LSB']['amp_corrected'];hc=ch2['carrier']['amp_corrected']/ch1['carrier']['amp_corrected'];hu=ch2['USB']['amp_corrected']/ch1['USB']['amp_corrected']
            row.update({'H_L':hl,'H_C':hc,'H_U':hu,'T_lower':hl/hc,'T_upper':hu/hc,'T_mean':(hl/hc+hu/hc)/2})
            prior=pd.concat([prior,pd.DataFrame([row])],ignore_index=True);prior.to_csv(REPEAT_CSV,index=False);done.add((fm,rep))
            say(f'{fm/1000:g}kHz rep{rep}: SNR CH1 L/U={ch1["LSB"]["snr"]:.1f}/{ch1["USB"]["snr"]:.1f}dB, CH2 L/U={ch2["LSB"]["snr"]:.1f}/{ch2["USB"]["snr"]:.1f}dB; T={hl/hc:.5f}/{hu/hc:.5f}; bins match={all(row[f"{c}_same_bin_CH1_CH2"] for c in ("LSB","carrier","USB"))}')
    if len(prior)!=len(FM_LIST)*REPS:raise RuntimeError(f'Sweep incomplete: {len(prior)}/{len(FM_LIST)*REPS} records')
    # Aggregate means and repeatability in a one-row-per-frequency delivery CSV.
    summaries=[]
    metric_cols=['H_L','H_C','H_U','T_lower','T_upper','T_mean']
    for fm,g in prior.groupby('fm_Hz'):
        row={'fm_Hz':int(fm),'repeats':len(g),'amplitude_programmed_Vpp':AMP,'carrier_programmed_Hz':FC,'depth_percent':DEPTH,'record_points':N,'sample_rate_Hz_mean':g.sample_rate_Hz.mean(),'record_duration_s_mean':g.record_duration_s.mean(),'fft_bin_Hz_mean':g.fft_bin_Hz.mean(),'all_captures_unclipped':not bool(g.CH1_clipped.any() or g.CH2_clipped.any()),'all_component_channel_bins_match':bool(g[[f'{c}_same_bin_CH1_CH2' for c in ('LSB','carrier','USB')]].all().all())}
        for c in metric_cols:
            row[c+'_mean']=g[c].mean();row[c+'_sd']=g[c].std(ddof=1);row[c+'_relative_sd_percent']=100*g[c].std(ddof=1)/abs(g[c].mean()) if g[c].mean()!=0 else math.nan
        for ch in (1,2):
            for comp in ('LSB','carrier','USB'):
                for metric in ('amp_corrected','noise_corrected','snr','frequency'):
                    col=f'CH{ch}_{comp}_{metric}'
                    row[col+'_mean']=g[col].mean();row[col+'_sd']=g[col].std(ddof=1)
        for ch in (1,2):
            for comp in ('LSB','USB'):
                col=f'CH{ch}_{comp}_snr';row[col+'_min']=g[col].min()
        summaries.append(row)
    summary=pd.DataFrame(summaries).sort_values('fm_Hz');summary.to_csv(SUMMARY_CSV,index=False)
    # SNR plot: both sidebands on both channels.
    fig,ax=plt.subplots(figsize=(9,5.8))
    for ch,comp,color,style in [(1,'LSB','#2563eb','-'),(1,'USB','#60a5fa','--'),(2,'LSB','#b91c1c','-'),(2,'USB','#f87171','--')]:
        y=f'CH{ch}_{comp}_snr_mean';sd=f'CH{ch}_{comp}_snr_sd';ax.errorbar(summary.fm_Hz/1000,summary[y],yerr=summary[sd],marker='o',ms=4,lw=1.3,capsize=2,color=color,ls=style,label=f'CH{ch} {comp}')
    ax.axhline(20,color='k',ls=':',lw=1,label='20 dB working adequacy threshold');ax.set(xlabel='Modulation frequency (kHz)',ylabel='Peak-bin / local median FFT floor (dB)',title='0.2 Vpp low-drive sideband SNR (mean ± 1 SD, n=3)');ax.grid(alpha=.25);ax.legend(ncol=2);fig.tight_layout();fig.savefig(ROOT/'low_drive_snr.png',dpi=180);plt.close(fig)
    # T mean +/- 1 SD and relative repeatability in two panels.
    fig,axs=plt.subplots(2,1,figsize=(9,7.2),sharex=True)
    x=summary.fm_Hz/1000
    for col,label,color in [('T_lower','T lower','#2563eb'),('T_upper','T upper','#b91c1c')]:axs[0].errorbar(x,summary[col+'_mean'],yerr=summary[col+'_sd'],marker='o',ms=4,lw=1.2,capsize=2,color=color,label=label)
    axs[0].set_ylabel('Normalized transfer T');axs[0].grid(alpha=.25);axs[0].legend();axs[0].set_xscale('log')
    axs[1].plot(x,summary.T_lower_relative_sd_percent,'o-',label='T lower relative SD',color='#2563eb');axs[1].plot(x,summary.T_upper_relative_sd_percent,'s-',label='T upper relative SD',color='#b91c1c');axs[1].set_xscale('log');axs[1].set(xlabel='Modulation frequency (kHz)',ylabel='Repeatability (relative SD, %)');axs[1].grid(alpha=.25);axs[1].legend();fig.suptitle('Three-capture sideband repeatability at 0.2 Vpp');fig.tight_layout();fig.savefig(ROOT/'low_drive_sideband_repeatability.png',dpi=180);plt.close(fig)
    say(f'COMPLETE: wrote {SUMMARY_CSV.name}; n={len(summary)} modulation frequencies x {REPS} repeats')
finally:
    if sc is not None:
        try:sc.write(':TRIG:RUN');sc.write(':TIM:SCAL 5E-2');sc.write(':CHAN1:SCAL 2E-1');sc.write(':CHAN2:SCAL 1E-1');sc.close()
        except Exception as e:say(f'Scope restore warning: {e!r}')
    if gn is not None:
        try:gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,FRQ,493900HZ,AMP,5V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ');say('FINAL generator '+gn.query('C2:BSWV?').strip()+'; '+gn.query('C2:MDWV?').strip()+'; '+gn.query('C2:OUTP?').strip());gn.close()
        except Exception as e:say(f'Generator restore warning: {e!r}')
    if rm is not None:rm.close()
