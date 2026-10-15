
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import time, pandas as pd, pyvisa
import characterize_rlc as c
S=_bench.scope_address;G=_bench.generator_address
rm=pyvisa.ResourceManager();sc=rm.open_resource(S);gn=rm.open_resource(G);sc.timeout=120000;gn.timeout=8000;sc.chunk_size=4*1024*1024
rows=[]
try:
 gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,FRQ,498400HZ,AMP,0.1V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ');sc.write(':TIM:SCAL 2E-4');time.sleep(.3)
 for name,s1,s2,n in [('5mV_5mV',.005,.005,5),('5mV_2mV',.005,.002,5),('2mV_2mV',.002,.002,3)]:
  sc.write(f':CHAN1:SCAL {s1}');sc.write(f':CHAN2:SCAL {s2}')
  for i in range(n):
   r=c.measure(sc,gn,498400.0,'range_repeat');r.update({'scale_condition':name,'repeat':i+1});rows.append(r)
  q=pd.DataFrame([x for x in rows if x['scale_condition']==name]);print(name,'Hmean',q.H_mag.mean(),'sd',q.H_mag.std(ddof=1),'Vin',q.vin_Vpk.mean(),'Vout',q.vout_Vpk.mean(),'phase',q.H_phase_deg.mean(),flush=True)
 pd.DataFrame(rows).to_csv('scope_scale_crosscheck_0p1vpp.csv',index=False)
finally:
 try:gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,FRQ,493900HZ,AMP,5V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ');sc.write(':CHAN1:SCAL 2E-1');sc.write(':CHAN2:SCAL 1E-1');sc.write(':TIM:SCAL 5E-2');sc.write(':TRIG:RUN')
 except Exception as e:print('restore',repr(e))
 sc.close();gn.close();rm.close()
