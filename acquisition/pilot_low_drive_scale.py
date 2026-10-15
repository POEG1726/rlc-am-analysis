
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import time, pyvisa
import characterize_rlc as c
rm=pyvisa.ResourceManager(); sc=rm.open_resource(c.SCOPE_ADDR);gn=rm.open_resource(c.GEN_ADDR);sc.timeout=120000;gn.timeout=8000;sc.chunk_size=4*1024*1024
try:
 gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,FRQ,497150HZ,AMP,0.1V,OFST,0V');gn.write('C2:OUTP ON,LOAD,HZ')
 sc.write(':TIM:SCAL 2E-4');sc.write(':CHAN1:SCAL 2E-3');sc.write(':CHAN2:SCAL 2E-3');time.sleep(.3)
 print('scales',sc.query(':CHAN1:SCAL?').strip(),sc.query(':CHAN2:SCAL?').strip())
 r=c.measure(sc,gn,497150.0,'pilot');print({k:r[k] for k in ['frequency_measured_Hz','vin_Vpk','vout_Vpk','H_mag','H_phase_deg','CH1_code_min','CH1_code_max','CH2_code_min','CH2_code_max','CH1_vdiv_V_div','CH2_vdiv_V_div']})
finally:
 gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,FRQ,493900HZ,AMP,5V,OFST,0V');sc.write(':TRIG:RUN');sc.write(':TIM:SCAL 5E-2');sc.write(':CHAN1:SCAL 2E-1');sc.write(':CHAN2:SCAL 1E-1');sc.close();gn.close();rm.close()
