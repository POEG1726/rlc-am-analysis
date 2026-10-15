
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import time
from pathlib import Path
import pandas as pd, pyvisa
import characterize_rlc as c
ROOT=_bench.workspace; rm=pyvisa.ResourceManager();sc=rm.open_resource(c.SCOPE_ADDR);gn=rm.open_resource(c.GEN_ADDR)
sc.timeout=120000;gn.timeout=8000;sc.chunk_size=4*1024*1024
rows=[]
try:
    gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,AMP,1.2V,OFST,0V,FRQ,500000HZ')
    sc.write(':TIM:SCAL 2E-4');time.sleep(.2)
    for f in [499000.0,500000.0,501000.0]:
        row=c.measure(sc,gn,f,'low_amplitude_pure_sine');rows.append(row)
        print(f'{f:.0f} Hz: Vin={row["vin_Vpk"]:.6f} Vpk Vout={row["vout_Vpk"]:.6f} Vpk H={row["H_mag"]:.7f} phase={row["H_phase_deg"]:+.4f} deg',flush=True)
    df=pd.DataFrame(rows);df.to_csv(ROOT/'amplitude_linearity_check.csv',index=False)
finally:
    try:
        gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,AMP,5V,OFST,0V,FRQ,493900HZ');sc.write(':TRIG:RUN');sc.write(':TIM:SCAL 5E-2');sc.close();gn.close();rm.close()
    except Exception:pass
