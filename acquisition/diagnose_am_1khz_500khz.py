
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json, re, struct, time
from pathlib import Path
import numpy as np, pandas as pd, pyvisa
ROOT=_bench.workspace
SCOPE=_bench.scope_address; GEN=_bench.generator_address
N=5_000_000; FM=1000.0; DEPTH=50.0; FC=500000.0
log=[]
def say(x): print(x,flush=True);log.append(str(x))
def block(inst,cmd):
    inst.write(cmd); raw=inst.read_raw(); p=raw.find(b'#')
    if p<0: raise RuntimeError(f'No block: {cmd}: {raw[:60]!r}')
    nd=int(chr(raw[p+1])); k=int(raw[p+2:p+2+nd]); data=raw[p+2+nd:p+2+nd+k]
    if len(data)!=k: raise RuntimeError(f'block truncated {cmd}: {len(data)}/{k}')
    return data
def capture(sc,ch):
    sc.write(f':WAV:SOUR C{ch}');sc.write(':WAV:WIDT WORD');sc.write(':WAV:BYT LSB');sc.write(':WAV:INT 1');sc.write(':WAV:STAR 0');sc.write(f':WAV:POIN {N}')
    pre=block(sc,':WAV:PRE?');vdiv,offset,cpd,dt=[struct.unpack_from('<f',pre,p)[0] for p in (156,160,164,176)]
    codes=np.frombuffer(block(sc,':WAV:DATA?'),dtype='<i2');v=codes.astype(float)*(vdiv/cpd)-offset
    return v,dt,dict(vdiv_V_div=vdiv,offset_V=offset,codes_per_div=cpd,n=len(codes),min_code=int(codes.min()),max_code=int(codes.max()))
def fft_peak(freq,amp,target,width=100):
    ids=np.flatnonzero(np.abs(freq-target)<=width); k=int(ids[np.argmax(amp[ids])]);delta=0
    if 0<k<len(amp)-1:
        a,b,c=np.log(np.maximum(amp[k-1:k+2],1e-300)); den=a-2*b+c
        if den:delta=.5*(a-c)/den
    return (k+delta)*(freq[1]-freq[0]),amp[k],k
try:
    rm=pyvisa.ResourceManager();sc=rm.open_resource(SCOPE);gn=rm.open_resource(GEN);sc.timeout=120000;gn.timeout=10000;sc.chunk_size=16*1024*1024
    say('Scope '+sc.query('*IDN?').strip());say('Generator '+gn.query('*IDN?').strip())
    sc.write(':TIM:SCAL 5E-2'); time.sleep(.3) # 50 ms total record at 5 ms/div
    gn.write('C2:BSWV WVTP,SINE,FRQ,500000HZ,AMP,5V,OFST,0V');time.sleep(.1)
    gn.write('C2:MDWV STATE,ON,AM,MDSP,SINE,SRC,INT,FRQ,1000HZ,DEPTH,50');time.sleep(.2)
    md=gn.query('C2:MDWV?').strip(); say('MDWV '+md)
    if 'STATE,ON' not in md.upper() or 'FRQ,1000HZ' not in md.upper(): raise RuntimeError('Generator AM configuration not confirmed')
    say('BSWV '+gn.query('C2:BSWV?').strip());say('OUTP '+gn.query('C2:OUTP?').strip())
    say(f'ACQ points={sc.query(":ACQ:POIN?").strip()} rate={sc.query(":ACQ:SRAT?").strip()} timebase={sc.query(":TIM:SCAL?").strip()}')
    say(f'CH impedance {sc.query(":CHAN1:IMP?").strip()} / {sc.query(":CHAN2:IMP?").strip()}')
    sc.write(':TRIG:RUN');time.sleep(.6);sc.write(':TRIG:STOP');time.sleep(.05)
    vin,dt1,meta1=capture(sc,1);vout,dt2,meta2=capture(sc,2)
    if len(vin)!=N or len(vout)!=N or not np.isclose(dt1,dt2):raise RuntimeError('Bad paired waveform record')
    fs=1/dt1; n=len(vin); duration=n/fs; df=fs/n
    say(f'Acquired N={n}, dt={dt1:.12g}, fs={fs:.6f}, duration={duration:.6f}s, df={df:.6f}Hz')
    np.savez_compressed(ROOT/'am_1khz_diagnostic_500khz.npz',vin=vin.astype(np.float32),vout=vout.astype(np.float32),dt_s=dt1,fc_set_Hz=FC,fm_set_Hz=FM)
    result={'frequency_component':'','Vin_Vpk':0,'Vout_Vpk':0,'H':0,'FFT_freq_in_Hz':0,'FFT_freq_out_Hz':0,'FFT_H':0}
    x1=vin-vin.mean();x2=vout-vout.mean();win=np.hanning(n);cg=win.mean();freq=np.fft.rfftfreq(n,dt1)
    a1=2*np.abs(np.fft.rfft(x1*win))/(n*cg);a2=2*np.abs(np.fft.rfft(x2*win))/(n*cg);a1[0]/=2;a2[0]/=2
    rows=[]
    for label,target in [('LSB',FC-FM),('carrier',FC),('USB',FC+FM)]:
        f1,v1,k1=fft_peak(freq,a1,target);f2,v2,k2=fft_peak(freq,a2,target)
        # Existing method selects independently nearest local maxima in each channel.
        rows.append({'component':label,'target_Hz':target,'Vin_fft_Vpk':v1,'Vin_fft_peak_Hz':f1,'Vout_fft_Vpk':v2,'Vout_fft_peak_Hz':f2,'H_fft':v2/v1,'in_bin':k1,'out_bin':k2})
        say(f'{label}: Vin={v1:.7f} at {f1:.3f} Hz; Vout={v2:.7f} at {f2:.3f} Hz; H={v2/v1:.7f}; bins={k1}/{k2}')
    h=np.array([r['H_fft'] for r in rows]); tl,tc,tu=h; say(f'FFT normalize Tlower={tl/tc:.7f}, Tupper={tu/tc:.7f}, Tmean={(tl+tu)/(2*tc):.7f}')
    dfout=pd.DataFrame(rows);dfout.to_csv(ROOT/'am_1khz_diagnostic_500khz.csv',index=False)
    # Compare old method to a common-bin ratio and local quadratic peak height from both channels.
    for r in rows:
        target=r['target_Hz'];idx=int(np.argmin(abs(freq-target)))
        say(f'{r["component"]} same-bin ratio at {freq[idx]:.1f}Hz = {a2[idx]/a1[idx]:.7f}')
    say('CH1 metadata '+json.dumps(meta1));say('CH2 metadata '+json.dumps(meta2))
    # AC amplitude envelope ranges and rail checks.
    say(f'CH1 raw range {vin.min():.6f}..{vin.max():.6f}V, CH2 {vout.min():.6f}..{vout.max():.6f}V')
    say(f'CH1 DC {vin.mean():.6g} VrmsAC {np.std(vin):.6g}; CH2 DC {vout.mean():.6g} VrmsAC {np.std(vout):.6g}')
    say('Wrote compressed paired native waveform vectors to am_1khz_diagnostic_500khz.npz')
    (ROOT/'am_1khz_diagnostic_500khz.log').write_text('\n'.join(log)+'\n',encoding='utf-8')
finally:
    try:
        if 'sc' in locals():sc.write(':TRIG:RUN');sc.write(':TIM:SCAL 5E-2');sc.close()
        if 'gn' in locals():
            gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV FRQ,493900HZ');gn.close()
        if 'rm' in locals():rm.close()
    except Exception:pass
