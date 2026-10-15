"""A/B/A scope-range and run-offset check at low-drive 1 kHz AM."""

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json,struct,time
from pathlib import Path
import numpy as np,pandas as pd,pyvisa
import characterize_rlc as base
ROOT=_bench.workspace;N=5_000_000;FC=498337.4;FM=1000;AMP=.2;DEPTH=50
SCOPE=base.SCOPE_ADDR;GEN=base.GEN_ADDR
cal=json.loads((ROOT/'vertical_scale_corrections.json').read_text(encoding='utf-8'))['correction_factors_to_reference']
def fac(ch,scale):return next(float(v) for k,v in cal[f'CH{ch}'].items() if np.isclose(float(k),scale,rtol=0,atol=1e-9))
def block(i,c):
 i.write(c);raw=i.read_raw();p=raw.find(b'#')
 if p<0:raise RuntimeError(f'No block: {c} {raw[:50]!r}')
 nd=int(chr(raw[p+1]));st=p+2+nd;n=int(raw[p+2:p+2+nd]);x=raw[st:st+n]
 if len(x)!=n:raise RuntimeError(f'truncated {c}: {len(x)}/{n}')
 return x
def cap(s,ch):
 s.write(f':WAV:SOUR C{ch}');s.write(':WAV:WIDT WORD');s.write(':WAV:BYT LSB');s.write(':WAV:INT 1');s.write(':WAV:STAR 0');s.write(f':WAV:POIN {N}')
 p=block(s,':WAV:PRE?');vd,off,cpd,dt=[struct.unpack_from('<f',p,x)[0] for x in(156,160,164,176)];c=np.frombuffer(block(s,':WAV:DATA?'),dtype='<i2');v=c.astype(float)*vd/cpd-off
 return v,dt,{'scale':vd,'min':int(c.min()),'max':int(c.max()),'clipped':bool(c.min()<=-32768 or c.max()>=32767)}
def lines(x,dt):
 n=len(x);w=np.hanning(n);a=2*np.abs(np.fft.rfft((x-x.mean())*w))/(n*w.mean());a[0]/=2;f=np.fft.rfftfreq(n,dt);res=[]
 for target in (FC-FM,FC,FC+FM):
  ids=np.flatnonzero(np.abs(f-target)<=100);k=int(ids[np.argmax(a[ids])]);res.append((int(k),float(f[k]),float(a[k])))
 return res
rm=pyvisa.ResourceManager();sc=rm.open_resource(SCOPE);gn=rm.open_resource(GEN);sc.timeout=120000;gn.timeout=10000;sc.chunk_size=16*1024*1024
rows=[];log=[]
def say(x):print(x,flush=True);log.append(x)
try:
 say('Scope '+sc.query('*IDN?').strip());say('Generator '+gn.query('*IDN?').strip())
 gn.write('C2:MDWV STATE,OFF');gn.write(f'C2:BSWV WVTP,SINE,FRQ,{FC}HZ,AMP,{AMP}V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ');gn.write(f'C2:MDWV STATE,ON,AM,MDSP,SINE,SRC,INT,FRQ,{FM}HZ,DEPTH,{DEPTH}')
 say('AM readback '+gn.query('C2:MDWV?').strip());say('Carrier readback '+gn.query('C2:BSWV?').strip())
 sc.write(':TIM:SCAL 5E-2')
 for stage,(s1,s2) in [('A1_original',(0.02,0.01)),('B_wider',(0.05,0.02)),('A2_return',(0.02,0.01))]:
  sc.write(f':CHAN1:SCAL {s1}');sc.write(f':CHAN2:SCAL {s2}');time.sleep(.5)
  say(f'{stage}: readback scales={sc.query(":CHAN1:SCAL?").strip()}/{sc.query(":CHAN2:SCAL?").strip()} V/div')
  for rep in range(1,4):
   sc.write(':TRIG:RUN');time.sleep(.65);sc.write(':TRIG:STOP');time.sleep(.05)
   x1,dt1,m1=cap(sc,1);x2,dt2,m2=cap(sc,2)
   if len(x1)!=N or len(x2)!=N or not np.isclose(dt1,dt2):raise RuntimeError('bad pair')
   l1=lines(x1,dt1);l2=lines(x2,dt2);corr1=fac(1,s1);corr2=fac(2,s2)
   hs=[]
   for i,(a,b) in enumerate(zip(l1,l2)):
    h=(b[2]*corr2)/(a[2]*corr1);hs.append(h)
   tl,tu=hs[0]/hs[1],hs[2]/hs[1];row={'stage':stage,'repeat':rep,'fm_Hz':FM,'carrier_programmed_Hz':FC,'amp_Vpp':AMP,'CH1_scale_V_div':s1,'CH2_scale_V_div':s2,'CH1_scale_correction':corr1,'CH2_scale_correction':corr2,'CH1_clipped':m1['clipped'],'CH2_clipped':m2['clipped'],'CH1_code_min':m1['min'],'CH1_code_max':m1['max'],'CH2_code_min':m2['min'],'CH2_code_max':m2['max'],'CH1_adc_fraction':max(abs(m1['min']),abs(m1['max']))/32768,'CH2_adc_fraction':max(abs(m2['min']),abs(m2['max']))/32768,'Vin_carrier_Vpk_corrected':l1[1][2]*corr1,'Vout_carrier_Vpk_corrected':l2[1][2]*corr2,'H_L':hs[0],'H_C':hs[1],'H_U':hs[2],'T_lower':tl,'T_upper':tu,'T_mean':(tl+tu)/2,'all_three_bins_match':all(a[0]==b[0] for a,b in zip(l1,l2)),'all_peak_frequencies_Hz':';'.join(f'{z[1]:.6f}' for z in l1)}
   rows.append(row);pd.DataFrame(rows).to_csv(ROOT/'low_drive_1khz_range_ab.csv',index=False)
   say(f'{stage} rep{rep}: codes CH1 {m1["min"]}..{m1["max"]} CH2 {m2["min"]}..{m2["max"]}, clipped={m1["clipped"]}/{m2["clipped"]}, T={tl:.6f}/{tu:.6f}, Hc={hs[1]:.6f}, bin match={row["all_three_bins_match"]}')
 pd.DataFrame(rows).to_csv(ROOT/'low_drive_1khz_range_ab.csv',index=False)
 summary=pd.DataFrame(rows).groupby('stage').agg({c:['mean','std','min','max'] for c in ['T_lower','T_upper','T_mean','H_C','Vin_carrier_Vpk_corrected','Vout_carrier_Vpk_corrected','CH1_adc_fraction','CH2_adc_fraction']})
 summary.to_csv(ROOT/'low_drive_1khz_range_ab_summary.csv');say('AB summary\n'+summary.to_string())
 (ROOT/'low_drive_1khz_range_ab.log').write_text('\n'.join(log)+'\n',encoding='utf-8')
finally:
 try:gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,FRQ,493900HZ,AMP,5V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ');say('RESTORE '+gn.query('C2:BSWV?').strip()+'; '+gn.query('C2:MDWV?').strip()+'; '+gn.query('C2:OUTP?').strip())
 except Exception as e:say(f'generator restore error {e!r}')
 try:sc.write(':TRIG:RUN');sc.write(':TIM:SCAL 5E-2');sc.write(':CHAN1:SCAL 2E-1');sc.write(':CHAN2:SCAL 1E-1');say('RESTORE scope '+sc.query(':TIM:SCAL?').strip()+'; '+sc.query(':CHAN1:SCAL?').strip()+'; '+sc.query(':CHAN2:SCAL?').strip())
 except Exception as e:say(f'scope restore error {e!r}')
 sc.close();gn.close();rm.close();(ROOT/'low_drive_1khz_range_ab.log').write_text('\n'.join(log)+'\n',encoding='utf-8')
