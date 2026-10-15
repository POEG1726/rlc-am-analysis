
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import time, json
from pathlib import Path
import numpy as np, pandas as pd, pyvisa
import characterize_rlc as c
ROOT=_bench.workspace
rm=pyvisa.ResourceManager(); sc=rm.open_resource(c.SCOPE_ADDR); gn=rm.open_resource(c.GEN_ADDR)
sc.timeout=120000; gn.timeout=8000; sc.chunk_size=4*1024*1024
rows=[]
try:
    md=gn.query('C2:MDWV?').strip()
    if 'STATE,OFF' not in md.upper(): raise RuntimeError('Refinement requires AM OFF')
    sc.write(':TIM:SCAL 2E-4')
    # Coarse result places peak at 493.5 kHz programmed; span 2 kHz, 100 Hz steps.
    for i,f in enumerate(np.arange(491500.0,495500.1,100.0)):
        row=c.measure(sc,gn,float(f),'refined'); rows.append(row)
        if i%10==0: print(f'{i+1}/41 f_set={f:.0f} H={row["H_mag"]:.7f} phase={row["H_phase_deg"]:+.4f}deg',flush=True)
    ref=pd.DataFrame(rows); ref.to_csv(ROOT/'resonance_refinement.csv',index=False)
    all_df=pd.concat([pd.read_csv(ROOT/'resonance_scan.csv'),ref],ignore_index=True)
    all_df.to_csv(ROOT/'resonance_scan.csv',index=False)
    # Peak and phase-zero estimates on commanded-frequency axis, plus polynomial refinement.
    top=ref.nlargest(5,'H_mag').sort_values('frequency_set_Hz')
    fset=ref.frequency_set_Hz.to_numpy(); hm=ref.H_mag.to_numpy(); idx=int(np.argmax(hm))
    p=np.polyfit(fset[idx-2:idx+3],hm[idx-2:idx+3],2); fpk=-p[1]/(2*p[0])
    phase=np.unwrap(np.deg2rad(ref.H_phase_deg.to_numpy()))
    iz=int(np.argmin(np.abs(ref.frequency_set_Hz-493900)))
    # use linear interpolation across nearest sign change in phase
    crossings=np.where(phase[:-1]*phase[1:]<=0)[0]
    izero=int(crossings[np.argmin(np.abs(ref.frequency_set_Hz.to_numpy()[crossings]-493900))])
    fphase=float(np.interp(0,[phase[izero],phase[izero+1]],[fset[izero],fset[izero+1]]))
    print('REFINED PEAK',fpk,'Hz programmed, Hmax',float(np.polyval(p,fpk)),'nearest point',float(fset[idx]),float(hm[idx]))
    print('PHASE ZERO programmed Hz',fphase)
    print('Top points\n',top[['frequency_set_Hz','frequency_measured_Hz','H_mag','H_phase_deg']].to_string(index=False))
    # Keep the measured grid peak as carrier set point, avoiding overclaiming sub-grid precision.
    target=float(fphase)
    gn.write(f'C2:BSWV FRQ,{target:.2f}HZ'); time.sleep(.15)
    print('GENERATOR FINAL',gn.query('C2:BSWV?').strip(),gn.query('C2:MDWV?').strip(),gn.query('C2:OUTP?').strip())
    # Update compact refinement diagnostics.
    old=json.loads((ROOT/'characterization_parameters.json').read_text(encoding='utf-8'))
    old.update({'f0_measured_peak_programmed_Hz':fpk,'Hmax_refined':float(np.polyval(p,fpk)),'f0_phase_zero_programmed_Hz':fphase,'carrier_recommended_programmed_Hz':target,
                'f0_peak_scope_frequency_Hz':float(np.interp(fpk,fset,ref.frequency_measured_Hz)), 'refinement_step_Hz':100.0,'refinement_n':len(ref)})
    (ROOT/'characterization_parameters.json').write_text(json.dumps(old,indent=2),encoding='utf-8')
    # Improve plot with finer points overlaid.
    fig_ax=None
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,6))
    broad=all_df[all_df.label.isin(['coarse','fine'])].sort_values('frequency_measured_Hz')
    ax.plot(broad.frequency_measured_Hz/1000,broad.H_mag,'o',ms=2,label='Broad measured sweep')
    ax.plot(ref.frequency_measured_Hz/1000,ref.H_mag,'-',lw=1.4,label='100 Hz peak refinement')
    ax.axvline(float(np.interp(fpk,fset,ref.frequency_measured_Hz))/1000,color='k',ls='--',alpha=.6,label=f'Peak {fpk/1000:.3f} kHz programmed')
    ax.set(xlabel='Measured frequency (kHz)',ylabel='|Vout/Vin|',title='Series RLC resistor transfer response'); ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'resonance_curve.png',dpi=180); plt.close(fig)
finally:
    try: sc.write(':TRIG:RUN'); sc.write(':TIM:SCAL 5E-2'); sc.close(); gn.close(); rm.close()
    except Exception: pass
