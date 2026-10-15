"""Five complete formal AM sweeps with non-overwriting raw waveform archives."""
from __future__ import annotations

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json, math, os, re, struct, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import pyvisa
import matplotlib.pyplot as plt

ROOT=_bench.workspace
BASE=ROOT/'formal_data'; RAW=BASE/'raw'; PROC=BASE/'processed'; FIG=BASE/'figures'; LOGDIR=BASE/'logs'
SCOPE_ADDR=_bench.scope_address
GEN_ADDR=_bench.generator_address
FM_LIST=[1000,2000,3000,5000,8000,10000,12000,15000,18000,20000,22000,25000,30000,35000,40000]
FC=498337.4; AMP=0.2; DEPTH=50.0; N=5_000_000; FS=50_000_000.0; DT=1/FS
SCALES={1:0.020,2:0.010}; TIMEBASE=0.010; SWEEPS=5; MIN_SNR_DB=20.0
for d in (RAW,PROC,FIG,LOGDIR): d.mkdir(parents=True,exist_ok=True)
for k in range(1,SWEEPS+1):(RAW/f'sweep_{k:02d}').mkdir(parents=True,exist_ok=True)
SUMMARY_MD=BASE/'formal_data_collection_summary.md'
if SUMMARY_MD.exists():raise FileExistsError(f'Completed collection exists; refusing to overwrite {SUMMARY_MD}')
if any((PROC/n).exists() for n in ('formal_am_raw.csv','formal_am_summary.csv','formal_am_end_checks.csv')):
    raise FileExistsError('Processed formal outputs already exist; refusing to overwrite.')

def now():return datetime.now(timezone.utc).isoformat(timespec='milliseconds')
def log(path,msg):
    line=f'[{now()}] {msg}';print(line,flush=True)
    with path.open('a',encoding='utf-8') as f:f.write(line+'\n')
def block(inst,cmd):
    inst.write(cmd);raw=inst.read_raw();p=raw.find(b'#')
    if p<0:raise RuntimeError(f'No IEEE block for {cmd}: {raw[:64]!r}')
    nd=int(chr(raw[p+1]));a=p+2;b=a+nd;n=int(raw[a:b].decode('ascii'));out=raw[b:b+n]
    if len(out)!=n:raise RuntimeError(f'Truncated {cmd}: {len(out)}/{n} bytes')
    return out
def read_channel(scope,ch):
    scope.write(f':WAV:SOUR C{ch}');pre=block(scope,':WAV:PRE?')
    if len(pre)<188:raise RuntimeError(f'CH{ch} waveform preamble too short')
    vdiv,offset,cpd,dt=[struct.unpack_from('<f',pre,p)[0] for p in (156,160,164,176)]
    codes=np.frombuffer(block(scope,':WAV:DATA?'),dtype='<i2').copy()
    m={'vdiv_V_div':vdiv,'offset_V':offset,'codes_per_div':cpd,'dt_s':dt,
       'points':len(codes),'code_min':int(codes.min()),'code_max':int(codes.max())}
    return codes,m
def volts(c,m):return c.astype(np.float64)*(m['vdiv_V_div']/m['codes_per_div'])-m['offset_V']
def line_info(freq,a,target):
    ids=np.flatnonzero(np.abs(freq-target)<=100.0)
    if not len(ids):raise RuntimeError(f'No FFT bin near {target:g} Hz')
    k=int(ids[np.argmax(a[ids])]);df=freq[1]-freq[0];delta=0.0
    if 0<k<len(a)-1:
        u,v,w=np.log(np.maximum(a[k-1:k+2],1e-300));den=u-2*v+w
        if den:delta=.5*(u-w)/den
    mask=(np.abs(freq-target)>=100)&(np.abs(freq-target)<=450);na=a[mask];nf=freq[mask]
    if len(na)<20:raise RuntimeError('Insufficient local FFT noise bins')
    med=float(np.median(na));outlier=na>med*10**(12/20);keep=np.ones(len(na),bool)
    for fsp in nf[outlier]:keep &= np.abs(nf-fsp)>20
    floor=float(np.median(na[keep])) if keep.sum()>=10 else med
    amp=float(a[k])
    return {'frequency_Hz':float(freq[k]+delta*df),'amplitude_Vpk':amp,'bin_index':k,
            'noise_floor_Vpk_per_bin':floor,'snr_dB':20*math.log10(max(amp,1e-30)/max(floor,1e-30)),
            'noise_bins_used':int(keep.sum()),'spur_bins_excluded':int((~keep).sum())}
def analyze(c1,c2,m1,m2,fm,dt):
    if len(c1)!=N or len(c2)!=N:raise RuntimeError(f'Expected {N} points, got {len(c1)}/{len(c2)}')
    if not np.isclose(dt,DT,rtol=0,atol=1e-15) or not np.isclose(m1['dt_s'],m2['dt_s'],rtol=0,atol=1e-15):
        raise RuntimeError(f'Sample interval mismatch: {dt}, {m1["dt_s"]}, {m2["dt_s"]}')
    n=len(c1);freq=np.fft.rfftfreq(n,dt);win=np.hanning(n);cg=win.mean()
    lines={}
    for ch,c,m in ((1,c1,m1),(2,c2,m2)):
        v=volts(c,m);a=2*np.abs(np.fft.rfft((v-v.mean())*win))/(n*cg);a[0]/=2
        for name,target in (('LSB',FC-fm),('carrier',FC),('USB',FC+fm)):
            lines[f'CH{ch}_{name}']=line_info(freq,a,target)
        del v,a
    bins={name:lines[f'CH1_{name}']['bin_index']==lines[f'CH2_{name}']['bin_index'] for name in ('LSB','carrier','USB')}
    h={name:lines[f'CH2_{name}']['amplitude_Vpk']/lines[f'CH1_{name}']['amplitude_Vpk'] for name in ('LSB','carrier','USB')}
    tl=h['LSB']/h['carrier'];tu=h['USB']/h['carrier']
    return {'lines':lines,'same_bin':bins,'H_L':h['LSB'],'H_C':h['carrier'],'H_U':h['USB'],
      'T_lower':tl,'T_upper':tu,'T_mean':(tl+tu)/2,
      'CH1_ADC_min':int(c1.min()),'CH1_ADC_max':int(c1.max()),'CH2_ADC_min':int(c2.min()),'CH2_ADC_max':int(c2.max()),
      'CH1_clipped':bool(c1.min()<=-32768 or c1.max()>=32767),'CH2_clipped':bool(c2.min()<=-32768 or c2.max()>=32767),
      'CH1_code_fraction':float(max(abs(c1.min()),abs(c1.max()))/32768),'CH2_code_fraction':float(max(abs(c2.min()),abs(c2.max()))/32768),
      'CH1_wave_min_V':float(volts(c1,m1).min()),'CH1_wave_max_V':float(volts(c1,m1).max()),
      'CH2_wave_min_V':float(volts(c2,m2).min()),'CH2_wave_max_V':float(volts(c2,m2).max())}
def save_unique(path,**data):
    if path.exists():raise FileExistsError(f'Never overwrite acquisition: {path}')
    tmp=path.with_name(path.stem+f'.partial_{os.getpid()}'+path.suffix)
    with tmp.open('xb') as f:np.savez_compressed(f,**data)
    if path.exists():tmp.unlink();raise FileExistsError(f'Existing acquisition preserved: {path}')
    os.replace(tmp,path)
def payload(c1,c2,m1,m2,dt,fm,sweep,typ,t0,t1,bswv,mdwv,cfg,a,accepted,reason=''):
    empty=np.array([],dtype=np.int16)
    return {'ch1_codes':empty if c1 is None else c1.astype(np.int16,copy=False),
      'ch2_codes':empty if c2 is None else c2.astype(np.int16,copy=False),
      'dt_s':np.asarray(np.nan if dt is None else dt),'sample_rate_Hz':np.asarray(np.nan if dt is None else 1/dt),
      'number_of_points':np.asarray(0 if c1 is None else len(c1)),
      'programmed_carrier_Hz':np.asarray(FC),'programmed_modulation_Hz':np.asarray(fm),
      'programmed_amplitude_Vpp':np.asarray(AMP),'AM_depth_percent':np.asarray(DEPTH),
      'CH1_vertical_scale_V_div':np.asarray(np.nan if m1 is None else m1['vdiv_V_div']),
      'CH2_vertical_scale_V_div':np.asarray(np.nan if m2 is None else m2['vdiv_V_div']),
      'CH1_voltage_offset_V':np.asarray(np.nan if m1 is None else m1['offset_V']),
      'CH2_voltage_offset_V':np.asarray(np.nan if m2 is None else m2['offset_V']),
      'CH1_codes_per_div':np.asarray(np.nan if m1 is None else m1['codes_per_div']),
      'CH2_codes_per_div':np.asarray(np.nan if m2 is None else m2['codes_per_div']),
      'timestamp_start_utc':np.asarray(t0),'timestamp_end_utc':np.asarray(t1),
      'sweep_number':np.asarray(sweep),'acquisition_type':np.asarray(typ),
      'generator_BSWV_readback':np.asarray(bswv),'generator_MDWV_readback':np.asarray(mdwv),
      'scope_configuration_readback':np.asarray(cfg),
      'CH1_code_min':np.asarray(-32768 if m1 is None else m1['code_min']),
      'CH1_code_max':np.asarray(32767 if m1 is None else m1['code_max']),
      'CH2_code_min':np.asarray(-32768 if m2 is None else m2['code_min']),
      'CH2_code_max':np.asarray(32767 if m2 is None else m2['code_max']),
      'accepted':np.asarray(bool(accepted)),'failure_reason':np.asarray(reason),
      'analysis_json':np.asarray(json.dumps(a or {},separators=(',',':')))}
def scope_config(s):
    return {'timebase_s_div':float(s.query(':TIM:SCAL?')),'sample_rate_Hz':float(s.query(':ACQ:SRAT?')),
      'points':int(float(s.query(':ACQ:POIN?'))),'CH1_scale_V_div':float(s.query(':CHAN1:SCAL?')),
      'CH2_scale_V_div':float(s.query(':CHAN2:SCAL?')),'CH1_impedance':s.query(':CHAN1:IMP?').strip(),
      'CH2_impedance':s.query(':CHAN2:IMP?').strip(),'CH1_coupling':s.query(':CHAN1:COUP?').strip(),
      'CH2_coupling':s.query(':CHAN2:COUP?').strip(),'wave_points':int(float(s.query(':WAV:POIN?'))),
      'wave_interval':int(float(s.query(':WAV:INT?')))}
def verify_gen(g,fm):
    b=g.query('C2:BSWV?').strip();m=g.query('C2:MDWV?').strip();o=g.query('C2:OUTP?').strip()
    for pat,val in [(r'WVTP,SINE',b),(r'FRQ,498337\.4HZ',b),(r'AMP,0\.2V',b),(r'OFST,0V',b),
                    (r'STATE,ON',m),(r'MDSP,SINE',m),(r'SRC,INT',m),(fr'FRQ,{fm}HZ',m),(r'DEPTH,50',m),(r'OFST,0V',m),(r'ON,LOAD,HZ',o)]:
        if not re.search(pat,val.upper()):raise RuntimeError(f'Readback mismatch for {pat}: {val}')
    return b,m,o
def check_scope(s,c):
    if c['points']!=N or not np.isclose(c['sample_rate_Hz'],FS,rtol=0,atol=1) or not np.isclose(c['timebase_s_div'],TIMEBASE,rtol=0,atol=1e-12):
        raise RuntimeError(f'Scope acquisition setting mismatch: {c}')
    if not np.isclose(c['CH1_scale_V_div'],SCALES[1]) or not np.isclose(c['CH2_scale_V_div'],SCALES[2]):raise RuntimeError(f'Scope vertical scale mismatch: {c}')
    if c['CH1_impedance'].upper()!='ONEMEG' or c['CH2_impedance'].upper()!='ONEMEG':raise RuntimeError(f'Scope input impedance mismatch: {c}')
    if c['wave_points']!=N or c['wave_interval']!=1:raise RuntimeError(f'Wave transfer mismatch: {c}')
    if c['CH1_coupling'].upper()!='DC' or c['CH2_coupling'].upper()!='DC':raise RuntimeError(f'Scope coupling mismatch: {c}')
def accept_analysis(a):
    if a['CH1_clipped'] or a['CH2_clipped']:return 'ADC clipping'
    if not all(a['same_bin'].values()):return f'CH1/CH2 FFT-bin mismatch {a["same_bin"]}'
    side=[a['lines'][f'CH{ch}_{s}']['snr_dB'] for ch in (1,2) for s in ('LSB','USB')]
    if min(side)<MIN_SNR_DB:return f'sideband SNR below {MIN_SNR_DB} dB ({min(side):.2f} dB)'
    return ''
def load_existing(path,sweep,fm,typ):
    with np.load(path,allow_pickle=False) as d:
        if int(d['sweep_number'])!=sweep or int(d['programmed_modulation_Hz'])!=fm or str(d['acquisition_type'])!=typ or not bool(d['accepted']):
            raise RuntimeError(f'Existing acquisition metadata mismatch; preserving file and stopping: {path}')
        a=json.loads(str(d['analysis_json']))
        return a,str(d['timestamp_start_utc']),str(d['timestamp_end_utc'])
def save_rejected(folder,stem,kw):
    i=1
    while True:
        p=folder/f'{stem}_rejected_{i:02d}.npz'
        if not p.exists():save_unique(p,**kw);return p
        i+=1
def make_row(a,sweep,fm,typ,t0,t1,path):
    r={'sweep_number':sweep,'fm_Hz':fm,'acquisition_type':typ,'timestamp_start_utc':t0,'timestamp_end_utc':t1,
       'raw_waveform_filename':str(path.relative_to(BASE)),'programmed_carrier_Hz':FC,'programmed_amplitude_Vpp':AMP,
       'AM_depth_percent':DEPTH,'record_points':N,'sample_rate_Hz':FS,'record_duration_s':N/FS,'fft_bin_Hz':FS/N,
       'CH1_vertical_scale_V_div':SCALES[1],'CH2_vertical_scale_V_div':SCALES[2]}
    for k,v in a.items():
        if k not in ('lines','same_bin'):r[k]=v
    for name,v in a['same_bin'].items():r[f'{name}_same_bin_CH1_CH2']=v
    for name,item in a['lines'].items():
        for k,v in item.items():r[f'{name}_{k}']=v
    return r

rm=sc=gen=None
try:
    rm=pyvisa.ResourceManager();sc=rm.open_resource(SCOPE_ADDR);gen=rm.open_resource(GEN_ADDR)
    sc.timeout=180000;gen.timeout=15000;sc.chunk_size=16*1024*1024
    prelog=LOGDIR/'preflight.log'
    log(prelog,'PRE-FLIGHT start; prior characterization/pilot files are preserved.')
    log(prelog,'Scope ID '+sc.query('*IDN?').strip()+'; generator ID '+gen.query('*IDN?').strip())
    # Configure once before formal sweep 1.
    gen.write('C2:MDWV STATE,OFF');gen.write(f'C2:BSWV WVTP,SINE,FRQ,{FC}HZ,AMP,{AMP}V,OFST,0V')
    gen.write('C2:OUTP ON,LOAD,HZ');gen.write('C2:MDWV STATE,ON,AM,MDSP,SINE,SRC,INT,FRQ,1000HZ,DEPTH,50')
    for ch in (1,2):
        sc.write(f':CHAN{ch}:IMP ONEMeg');sc.write(f':CHAN{ch}:COUP DC');sc.write(f':CHAN{ch}:PROB 1');sc.write(f':CHAN{ch}:SCAL {SCALES[ch]:g}')
    sc.write(':ACQ:MDEP 5M');sc.write(f':TIM:SCAL {TIMEBASE:g}')
    sc.write(':WAV:WIDT WORD');sc.write(':WAV:BYT LSB');sc.write(':WAV:INT 1');sc.write(':WAV:STAR 0');sc.write(f':WAV:POIN {N}')
    time.sleep(.3);cfg=scope_config(sc);check_scope(sc,cfg);b,m,o=verify_gen(gen,1000)
    log(prelog,'Fixed preflight configuration '+json.dumps(cfg,sort_keys=True))
    log(prelog,'Generator readbacks BSWV='+b+'; MDWV='+m+'; OUTP='+o)
    # Verify the frequency-only command and return to the required first point.
    gen.write('C2:MDWV FRQ,2000HZ');time.sleep(.05);verify_gen(gen,2000)
    gen.write('C2:MDWV FRQ,1000HZ');time.sleep(.05);b,m,o=verify_gen(gen,1000)
    log(prelog,'Confirmed frequency-only MDWV update; returned to 1 kHz before preflight capture.')
    t0=now();sc.write(':TRIG:RUN');time.sleep(.30);sc.write(':TRIG:STOP');time.sleep(.025)
    c1,m1=read_channel(sc,1);c2,m2=read_channel(sc,2);t1=now();a=analyze(c1,c2,m1,m2,1000,m1['dt_s']);reason=accept_analysis(a)
    prevalid=(not reason and len(c1)==N and len(c2)==N)
    p=RAW/f'preflight_1kHz_{datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")}.npz'
    save_unique(p,**payload(c1,c2,m1,m2,m1['dt_s'],1000,0,'preflight',t0,t1,b,m,json.dumps(cfg,sort_keys=True),a,prevalid,reason))
    mins=[a['lines'][f'CH{ch}_{ss}']['snr_dB'] for ch in (1,2) for ss in ('LSB','USB')]
    log(prelog,f'Expected 1 kHz AM triplet preflight accepted={prevalid}; T={a["T_lower"]:.6f}/{a["T_upper"]:.6f}/{a["T_mean"]:.6f}; min sideband SNR={min(mins):.2f} dB; bins={a["same_bin"]}; clip={a["CH1_clipped"]}/{a["CH2_clipped"]}; file={p.relative_to(BASE)}')
    if not prevalid:raise RuntimeError(f'Preflight failed ({reason}); no formal sweep was started.')
    log(prelog,'PRE-FLIGHT passed. Beginning five sweeps; formal-phase instrument setting changes are modulation frequency only.')

    # Collect five independent, ordered sweeps and one 1 kHz end check after each 40 kHz point.
    for sweep in range(1,SWEEPS+1):
        folder=RAW/f'sweep_{sweep:02d}';lp=LOGDIR/f'sweep_{sweep:02d}.log'
        log(lp,f'SWEEP {sweep:02d} START')
        b,m,o=verify_gen(gen,1000);cfg=scope_config(sc);check_scope(sc,cfg)
        log(lp,'Instrument IDs scope='+sc.query('*IDN?').strip()+'; generator='+gen.query('*IDN?').strip())
        log(lp,'Start settings BSWV='+b+'; MDWV='+m+'; OUTP='+o+'; SCOPE='+json.dumps(cfg,sort_keys=True))
        plan=[(fm,'formal',f'fm_{fm//1000:02d}kHz') for fm in FM_LIST]+[(1000,'endcheck','endcheck_01kHz')]
        for fm,typ,stem in plan:
            # Only change the modulation frequency during a sweep.
            gen.write(f'C2:MDWV FRQ,{fm}HZ');time.sleep(.12);b,m,o=verify_gen(gen,fm)
            log(lp,f'fm set/readback {fm} Hz ({typ}); MDWV={m}')
            path=folder/(stem+'.npz')
            if path.exists():
                a,t0,t1=load_existing(path,sweep,fm,typ)
                log(lp,f'RESUME existing accepted record without overwrite {path.name}: T={a["T_lower"]:.6f}/{a["T_upper"]:.6f}/{a["T_mean"]:.6f}')
                continue
            ok=False
            for attempt in range(1,4):
                t0=now();c1=c2=m1=m2=None;aa=None;dt=None;t1='';reason=''
                try:
                    sc.write(':TRIG:RUN');time.sleep(.25);sc.write(':TRIG:STOP');time.sleep(.025)
                    c1,m1=read_channel(sc,1);c2,m2=read_channel(sc,2);t1=now();dt=m1['dt_s']
                    cfg_now=scope_config(sc);check_scope(sc,cfg_now)
                    aa=analyze(c1,c2,m1,m2,fm,dt);reason=accept_analysis(aa)
                    if reason:raise RuntimeError(reason)
                    ok=True
                except Exception as exc:
                    reason=repr(exc);t1=t1 or now()
                    rejected=save_rejected(folder,stem,payload(c1,c2,m1,m2,dt,fm,sweep,'rejected_'+typ,t0,t1,b,m,json.dumps(cfg,sort_keys=True),aa,False,reason))
                    log(lp,f'REJECTED {fm} Hz attempt {attempt}: {reason}; preserved {rejected.name}; reacquiring')
                    if 'Scope acquisition' in reason or 'Wave transfer' in reason or 'Scope vertical' in reason:raise
                    continue
                if ok:
                    save_unique(path,**payload(c1,c2,m1,m2,dt,fm,sweep,typ,t0,t1,b,m,json.dumps(cfg_now,sort_keys=True),aa,True))
                    sn=[aa['lines'][f'CH{ch}_{side}']['snr_dB'] for ch in (1,2) for side in ('LSB','USB')]
                    log(lp,f'VALID {fm} Hz attempt {attempt}: {path.name}; min SNR={min(sn):.2f} dB; bins={aa["same_bin"]}; ADC fractions CH1/CH2={aa["CH1_code_fraction"]:.3f}/{aa["CH2_code_fraction"]:.3f}; T={aa["T_lower"]:.6f}/{aa["T_upper"]:.6f}/{aa["T_mean"]:.6f}')
                    break
            if not ok:raise RuntimeError(f'No valid acquisition after 3 attempts for sweep {sweep}, fm={fm}')
        end=now();log(lp,f'SWEEP {sweep:02d} COMPLETE {end}; final end-check logged above.')

    # Build tables only from accepted, separately stored native-waveform records.
    formal_rows=[];end_rows=[];all_snr=[];rejected_files=[];all_clip=[];sweep_times=[]
    for sweep in range(1,SWEEPS+1):
        folder=RAW/f'sweep_{sweep:02d}'
        for fm in FM_LIST:
            path=folder/f'fm_{fm//1000:02d}kHz.npz'
            a,t0,t1=load_existing(path,sweep,fm,'formal')
            formal_rows.append(make_row(a,sweep,fm,'formal',t0,t1,path))
            for ch in (1,2):
                all_clip.append(a[f'CH{ch}_clipped'])
                for side in ('LSB','USB'):all_snr.append(a['lines'][f'CH{ch}_{side}']['snr_dB'])
        p0=folder/'fm_01kHz.npz';pe=folder/'endcheck_01kHz.npz'
        a0,t00,t01=load_existing(p0,sweep,1000,'formal');ae,te0,te1=load_existing(pe,sweep,1000,'endcheck')
        rd={'sweep_number':sweep,'initial_timestamp_start_utc':t00,'endcheck_timestamp_start_utc':te0,
            'initial_raw_waveform_filename':str(p0.relative_to(BASE)),'endcheck_raw_waveform_filename':str(pe.relative_to(BASE))}
        for key in ('T_lower','T_upper','T_mean'):
            x=float(a0[key]);y=float(ae[key]);rd[f'initial_{key}']=x;rd[f'endcheck_{key}']=y
            rd[f'{key}_drift']=y-x;rd[f'{key}_drift_percent']=100*(y-x)/x if x else np.nan
        end_rows.append(rd)
        files=list(folder.glob('*_rejected_*.npz'));rejected_files.extend(files)
        lp=LOGDIR/f'sweep_{sweep:02d}.log';log(lp,f'End-check result: initial T={a0["T_lower"]:.6f}/{a0["T_upper"]:.6f}/{a0["T_mean"]:.6f}; end T={ae["T_lower"]:.6f}/{ae["T_upper"]:.6f}/{ae["T_mean"]:.6f}; drift={rd["T_lower_drift"]:+.6g}/{rd["T_upper_drift"]:+.6g}/{rd["T_mean_drift"]:+.6g}.')
    if len(formal_rows)!=75 or len(end_rows)!=5:raise RuntimeError(f'Incomplete dataset: formal={len(formal_rows)}, endchecks={len(end_rows)}')
    raw_df=pd.DataFrame(formal_rows).sort_values(['sweep_number','fm_Hz'])
    end_df=pd.DataFrame(end_rows).sort_values('sweep_number')
    raw_out=PROC/'formal_am_raw.csv';sum_out=PROC/'formal_am_summary.csv';end_out=PROC/'formal_am_end_checks.csv'
    for p in (raw_out,sum_out,end_out):
        if p.exists():raise FileExistsError(f'Refusing to overwrite {p}')
    raw_df.to_csv(raw_out,index=False);end_df.to_csv(end_out,index=False)

    summary=[]
    for fm,g in raw_df.groupby('fm_Hz',sort=True):
        r={'fm_Hz':int(fm),'n_sweeps':int(len(g))}
        for key in ('T_lower','T_upper','T_mean'):
            vals=g[key].astype(float);mu=float(vals.mean());sd=float(vals.std(ddof=1))
            r[f'{key}_mean']=mu;r[f'{key}_sample_SD']=sd;r[f'{key}_relative_SD_percent']=100*sd/abs(mu) if mu else np.nan
        for key in ('H_L','H_C','H_U'):
            r[f'{key}_mean']=float(g[key].mean());r[f'{key}_sample_SD']=float(g[key].std(ddof=1))
        for col in ('CH1_LSB_snr_dB','CH1_USB_snr_dB','CH2_LSB_snr_dB','CH2_USB_snr_dB'):
            r[f'{col}_min']=float(g[col].min());r[f'{col}_mean']=float(g[col].mean())
        summary.append(r)
    sum_df=pd.DataFrame(summary);sum_df.to_csv(sum_out,index=False)

    # Figures: show individual sweeps alongside across-sweep means and sample SD.
    colors=plt.cm.viridis(np.linspace(.08,.92,SWEEPS))
    fig,ax=plt.subplots(figsize=(9,5.4))
    for sw,col in zip(range(1,SWEEPS+1),colors):
        g=raw_df[raw_df.sweep_number==sw];ax.plot(g.fm_Hz/1000,g.T_lower,'o-',color=col,label=f'Sweep {sw} lower',alpha=.8)
        ax.plot(g.fm_Hz/1000,g.T_upper,'s--',color=col,label=f'Sweep {sw} upper',alpha=.8)
    ax.set(xlabel='Modulation frequency (kHz)',ylabel='Normalized sideband transfer T',title='Lower and upper sideband transfer');ax.grid(True,alpha=.25);ax.legend(ncol=2,fontsize=7);fig.tight_layout();fig.savefig(FIG/'T_lower_upper.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.4));ax.errorbar(sum_df.fm_Hz/1000,sum_df.T_mean_mean,yerr=sum_df.T_mean_sample_SD,fmt='o-',capsize=3)
    ax.set(xlabel='Modulation frequency (kHz)',ylabel='T_mean',title='Mean normalized sideband transfer across five sweeps (±1 sample SD)');ax.grid(True,alpha=.25);fig.tight_layout();fig.savefig(FIG/'T_mean.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.4))
    for sw,col in zip(range(1,SWEEPS+1),colors):
        g=raw_df[raw_df.sweep_number==sw];ax.plot(g.fm_Hz/1000,g.T_mean,'o-',color=col,label=f'Sweep {sw}')
    ax.set(xlabel='Modulation frequency (kHz)',ylabel='T_mean',title='Sweep comparison');ax.grid(True,alpha=.25);ax.legend();fig.tight_layout();fig.savefig(FIG/'sweep_comparison.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.4))
    for key,label in (('T_lower','Lower'),('T_upper','Upper'),('T_mean','Mean')):
        # The CSV field is already expressed as percent (100 * SD / |mean|).
        ax.plot(sum_df.fm_Hz/1000,sum_df[f'{key}_relative_SD_percent'],'o-',label=label)
    ax.set(xlabel='Modulation frequency (kHz)',ylabel='Relative sample SD (%)',title='Across-sweep repeatability');ax.grid(True,alpha=.25);ax.legend();fig.tight_layout();fig.savefig(FIG/'repeatability.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.4))
    for key,label in (('T_lower','Lower'),('T_upper','Upper'),('T_mean','Mean')):
        ax.plot(end_df.sweep_number,end_df[f'{key}_drift_percent'],'o-',label=label)
    ax.axhline(0,color='black',lw=.8);ax.set(xlabel='Sweep',ylabel='End-check drift from initial 1 kHz (%)',title='Within-sweep 1 kHz drift');ax.grid(True,alpha=.25);ax.legend();fig.tight_layout();fig.savefig(FIG/'endcheck_drift.png',dpi=180);plt.close(fig)

    tmean_sd=float(sum_df.T_mean_sample_SD.max());worst_fm=int(sum_df.loc[sum_df.T_mean_sample_SD.idxmax(),'fm_Hz'])
    min_snr=float(min(all_snr));max_snr=float(max(all_snr));max_drift=float(end_df[[f'{k}_drift_percent' for k in ('T_lower','T_upper','T_mean')]].abs().to_numpy().max())
    # Flag anomalies using a transparent rule: a point is >3 pooled-SD from its frequency mean.
    anomalies=[]
    for _,r in raw_df.iterrows():
        sr=sum_df[sum_df.fm_Hz==r.fm_Hz].iloc[0];sd=float(sr.T_mean_sample_SD);mu=float(sr.T_mean_mean)
        if sd>0 and abs(float(r.T_mean)-mu)>3*sd:anomalies.append(f'sweep {int(r.sweep_number)} at {int(r.fm_Hz/1000)} kHz')
    accepted_count=len(raw_df)+len(end_df);rejected_count=len(rejected_files)
    starts=[x['timestamp_start_utc'] for x in formal_rows];ends=[x['timestamp_end_utc'] for x in formal_rows]
    md=[]
    md += ['# Formal AM data collection summary','',f'- Collection completed: {now()}','- Formal acquisitions: **75 / 75** (five sweeps × 15 frequencies).',f'- End checks: **{len(end_df)} / 5**; accepted acquisitions total: **{accepted_count} / 80**.',f'- Rejected attempts preserved: **{rejected_count}**.' if rejected_count else '- Rejected attempts preserved: **0**.',f'- Formal acquisition time span: {min(starts)} to {max(ends)}.','', '## Acquisition validity','',f'- Sideband SNR across CH1/CH2 and LSB/USB: **{min_snr:.1f} to {max_snr:.1f} dB** (local per-bin FFT estimate).',f'- ADC clipping in accepted records: **{"none" if not any(all_clip) else "detected"}**.', '- FFT-bin consistency: all accepted measurements passed the CH1/CH2 same-bin check; waveform data are native 5,000,000-point, 50 MSa/s captures.', '- No interpolation or acquisition-time averaging was applied. The summary uses the five separate sweeps.','', '## Agreement and repeatability','',f'- Largest sample SD of T_mean across modulation frequencies: **{tmean_sd:.6g}** at **{worst_fm/1000:g} kHz**.',f'- Largest absolute initial-to-end 1 kHz drift among T_lower, T_upper, and T_mean: **{max_drift:.4g}%**. See `processed/formal_am_end_checks.csv`.','', '| fm (kHz) | mean T_lower ± SD | mean T_upper ± SD | mean T_mean ± SD | max relative SD T_mean (%) |','|---:|---:|---:|---:|---:|']
    for _,r in sum_df.iterrows():md.append(f'| {r.fm_Hz/1000:g} | {r.T_lower_mean:.6f} ± {r.T_lower_sample_SD:.2g} | {r.T_upper_mean:.6f} ± {r.T_upper_sample_SD:.2g} | {r.T_mean_mean:.6f} ± {r.T_mean_sample_SD:.2g} | {r.T_mean_relative_SD_percent:.4g} |')
    md += ['', '## 1 kHz end checks','', '| Sweep | T_lower drift | T_upper drift | T_mean drift |','|---:|---:|---:|---:|']
    for _,r in end_df.iterrows():md.append(f'| {int(r.sweep_number)} | {r.T_lower_drift:+.6g} ({r.T_lower_drift_percent:+.4g}%) | {r.T_upper_drift:+.6g} ({r.T_upper_drift_percent:+.4g}%) | {r.T_mean_drift:+.6g} ({r.T_mean_drift_percent:+.4g}%) |')
    md += ['', '## Anomalies and suitability','',('- No point exceeded three times its across-sweep sample SD from that frequency mean.' if not anomalies else '- Points beyond three times their frequency sample SD: '+', '.join(anomalies)+'.'),'- Suitability: **suitable for IA analysis** if all five sweeps and end checks above are complete and no clipping was detected. Interpret the reported SDs as between-sweep repeatability; sideband SNR is estimated locally from FFT bins.','', '## Outputs','', '- `processed/formal_am_raw.csv` — one row per formal acquisition.','- `processed/formal_am_summary.csv` — per-frequency means and sample SDs across five sweeps.','- `processed/formal_am_end_checks.csv` — within-sweep initial-to-end drift.','- `raw/sweep_01` through `raw/sweep_05` — accepted native waveform NPZ files and any preserved rejects.','- `figures/` — requested plots.','']
    if any(all_clip):md[-5]='- Suitability: **review before IA analysis** because clipping was detected in accepted records.'
    SUMMARY_MD.write_text('\n'.join(md),encoding='utf-8')
    for sweep in range(1,SWEEPS+1):log(LOGDIR/f'sweep_{sweep:02d}.log','Post-collection processing completed; processed tables and requested figures saved.')
    log(prelog,f'ALL SWEEPS COMPLETE: {accepted_count} accepted acquisitions, {rejected_count} rejected attempts preserved; summary at {SUMMARY_MD}.')
finally:
    if sc is not None:
        try:sc.write(':TRIG:RUN')
        except Exception:pass
    if gen is not None:
        try:gen.close()
        except Exception:pass
    if sc is not None:
        try:sc.close()
        except Exception:pass
    if rm is not None:
        try:rm.close()
        except Exception:pass
