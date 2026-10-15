
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import time, pandas as pd, numpy as np, pyvisa
import characterize_rlc as c
S=_bench.scope_address;G=_bench.generator_address
CH1_PAIRS=[(.15,498358,.005,.010),(.4,498118,.010,.020),(.75,497757,.020,.050),(2.0,496483,.050,.100),(5.0,493894,.100,.200)]
CH2_PAIRS=[(.1,498378,.002,.005),(.3,498222,.005,.010),(.6,497917,.010,.020),(1.0,497482,.020,.050),(3.0,495598,.050,.100),(5.0,493894,.100,.200)]
rm=pyvisa.ResourceManager();sc=rm.open_resource(S);gn=rm.open_resource(G);sc.timeout=120000;gn.timeout=8000;sc.chunk_size=4*1024*1024
rows=[]
try:
 gn.write('C2:MDWV STATE,OFF');sc.write(':TIM:SCAL 2E-4')
 for channel,pairs in [(1,CH1_PAIRS),(2,CH2_PAIRS)]:
  for amp,f0,low,high in pairs:
   gn.write(f'C2:BSWV WVTP,SINE,FRQ,{f0:.1f}HZ,AMP,{amp:.6g}V,OFST,0V');time.sleep(.3)
   other_ch=2 if channel==1 else 1
   fixed={1:{.15:.005,.4:.02,.75:.02,2.0:.1,5.0:.2},2:{.1:.005,.3:.01,.6:.02,1.0:.05,3.0:.1,5.0:.1}}[channel][amp]
   for scale in [low,high]:
    scopes={channel:scale,other_ch:fixed};sc.write(f':CHAN1:SCAL {scopes[1]}');sc.write(f':CHAN2:SCAL {scopes[2]}');time.sleep(.2)
    # Discard one record after the vertical range change, then retain three repeats.
    c.measure(sc,gn,f0,'scale_settling')
    for rep in range(3):
     r=c.measure(sc,gn,f0,'scale_calibration');r.update({'tested_channel':channel,'amplitude_programmed_Vpp':amp,'frequency_set_Hz':f0,'tested_scale_V_div':scale,'fixed_other_channel_scale_V_div':fixed,'repeat':rep+1});rows.append(r)
   subset=pd.DataFrame([r for r in rows if r['tested_channel']==channel and r['amplitude_programmed_Vpp']==amp])
   means=subset.groupby('tested_scale_V_div').agg({'vin_Vpk':'mean','vout_Vpk':'mean','H_mag':'mean','H_phase_deg':'mean'})
   print(f'CH{channel} at {amp:g}Vpp range pair {low:g}/{high:g}:\n{means.to_string()}\n',flush=True)
 d=pd.DataFrame(rows);d.to_csv('scope_vertical_range_calibration.csv',index=False)
finally:
 try:gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,FRQ,493900HZ,AMP,5V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ');sc.write(':CHAN1:SCAL 2E-1');sc.write(':CHAN2:SCAL 1E-1');sc.write(':TIM:SCAL 5E-2');sc.write(':TRIG:RUN')
 except Exception as e:print('restore',repr(e))
 sc.close();gn.close();rm.close()
