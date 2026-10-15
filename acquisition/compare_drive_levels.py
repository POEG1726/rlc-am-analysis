
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import time
from pathlib import Path
import numpy as np, pandas as pd, pyvisa
import characterize_rlc as c
ROOT=_bench.workspace; rm=pyvisa.ResourceManager();sc=rm.open_resource(c.SCOPE_ADDR);gn=rm.open_resource(c.GEN_ADDR)
sc.timeout=120000;gn.timeout=8000;sc.chunk_size=4*1024*1024;rows=[]
try:
    gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,AMP,1.2V,OFST,0V,FRQ,500000HZ')
    sc.write(':TIM:SCAL 2E-4');time.sleep(.2)
    for f in np.arange(488000,506000.1,1000):
        r=c.measure(sc,gn,float(f),'low_drive_coarse');rows.append(r)
    d=pd.DataFrame(rows); d=d.sort_values('frequency_set_Hz')
    ph=np.unwrap(np.deg2rad(d.H_phase_deg.to_numpy())); fs=d.frequency_set_Hz.to_numpy()
    changes=np.where(ph[:-1]*ph[1:]<=0)[0]
    if not len(changes):raise RuntimeError('No phase zero crossing in low-drive coarse scan')
    j=int(changes[np.argmin(np.abs(fs[changes]-497000))]); fzero=float(np.interp(0,[ph[j],ph[j+1]],[fs[j],fs[j+1]]))
    print(f'LOW-DRIVE coarse phase-zero ~{fzero:.1f} Hz; max H at {d.loc[d.H_mag.idxmax(),"frequency_set_Hz"]:.0f} Hz')
    for f in np.arange(fzero-1500,fzero+1500.1,100):
        r=c.measure(sc,gn,float(f),'low_drive_fine');rows.append(r)
    out=pd.DataFrame(rows);out.to_csv(ROOT/'amplitude_comparison.csv',index=False)
    for label,g in out.groupby('label'):
        print(label,'n=',len(g),'Hrange=',g.H_mag.min(),g.H_mag.max(),'Vinpkrange=',g.vin_Vpk.min(),g.vin_Vpk.max())
    # phase-zero refined, peak magnitude, and table of comparable frequencies
    fine=out[out.label=='low_drive_fine'].sort_values('frequency_set_Hz');p= np.polyfit(fine.frequency_set_Hz,fine.H_phase_deg,1); z=-p[1]/p[0]
    pk=fine.loc[fine.H_mag.idxmax()]
    print(f'LOW-DRIVE phase zero fit {z:.2f} Hz; max measured H={pk.H_mag:.7f} at set {pk.frequency_set_Hz:.0f}Hz; phase={pk.H_phase_deg:.4f}deg')
    print(out[out.label=='low_drive_coarse'][['frequency_set_Hz','vin_Vpk','H_mag','H_phase_deg']].to_string(index=False))
    print('FINAL TARGET setup will restore 5 Vpp, AM off, 493900 Hz')
finally:
    try:
        gn.write('C2:MDWV STATE,OFF');gn.write('C2:BSWV WVTP,SINE,AMP,5V,OFST,0V,FRQ,493900HZ');sc.write(':TRIG:RUN');sc.write(':TIM:SCAL 5E-2');sc.close();gn.close();rm.close()
    except Exception:pass
