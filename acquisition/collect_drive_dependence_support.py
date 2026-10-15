"""Fresh AM-OFF resonance repeats at six drive amplitudes; no prior measurements are loaded."""
from __future__ import annotations

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json, math, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import pyvisa
import matplotlib.pyplot as plt
import characterize_drive_amplitude as method

ROOT=_bench.workspace
OUT=ROOT/'drive_dependence_support_20260930_run03'
AMPS=[0.2,1.0,2.0,3.0,4.0,5.0]
REPEATS=3
R_EXT=5.1
FS=50_000_000.0
N=100_000
SUPPORTED=method.SUPPORTED_SCALES
NOMINAL_F0=1/(2*np.pi*np.sqrt(22e-6*4.7e-9))
CSV=OUT/'drive_dependence_support.csv'
RAW=OUT/'drive_dependence_frequency_scans.csv'
HARM=OUT/'drive_dependence_harmonics.csv'
LOG=OUT/'drive_dependence_support.log'
SUMMARY=OUT/'drive_dependence_support_summary.md'
RUNLOG=[]
RANGE_CHECKS=[]

def log(s):
    line=f'[{datetime.now(timezone.utc).isoformat(timespec="milliseconds")}] {s}'
    print(line,flush=True);RUNLOG.append(line)
    if OUT.exists():
        with LOG.open('a',encoding='utf-8') as f:f.write(line+'\n')

def step_scale(scale,direction):
    if direction>0:
        a=SUPPORTED[SUPPORTED>scale*(1+1e-9)]
        return float(a[0]) if len(a) else float(scale)
    a=SUPPORTED[SUPPORTED<scale*(1-1e-9)]
    return float(a[-1]) if len(a) else float(scale)

def fresh_auto_range_at(scope,gen,amp_vpp,fcenter):
    """Find vertical ranges using only freshly captured waveforms in this run."""
    scales={1:0.2,2:0.1}
    for ch in (1,2):scope.write(f':CHAN{ch}:SCAL {scales[ch]:.9g}')
    row=None
    for attempt in range(1,11):
        row=method.base.measure(scope,gen,float(fcenter),'fresh_range_check')
        divs={1:max(abs(row['CH1_code_min']),abs(row['CH1_code_max']))/7680.0,
              2:max(abs(row['CH2_code_min']),abs(row['CH2_code_max']))/7680.0}
        clipped={1:row['CH1_code_min']<=-32768 or row['CH1_code_max']>=32767,
                 2:row['CH2_code_min']<=-32768 or row['CH2_code_max']>=32767}
        changed=False
        for ch in (1,2):
            if clipped[ch] or divs[ch]>3.35:
                new=step_scale(scales[ch],1)
            elif divs[ch]<1.15:
                new=step_scale(scales[ch],-1)
            else:
                new=scales[ch]
            if new!=scales[ch]:
                scales[ch]=new;scope.write(f':CHAN{ch}:SCAL {new:.9g}');changed=True
        if not changed:break
    scales={ch:float(scope.query(f':CHAN{ch}:SCAL?')) for ch in (1,2)}
    divs={1:max(abs(row['CH1_code_min']),abs(row['CH1_code_max']))/7680.0,
          2:max(abs(row['CH2_code_min']),abs(row['CH2_code_max']))/7680.0}
    if any(divs[ch]>3.55 for ch in (1,2)):
        raise RuntimeError(f'Fresh auto-range remained too close to display limit: divisions={divs}, scales={scales}')
    RANGE_CHECKS.append({'amplitude_programmed_Vpp':amp_vpp,'frequency_set_Hz':fcenter,
                         'CH1_scale_V_div':scales[1],'CH2_scale_V_div':scales[2],
                         'CH1_peak_divisions':divs[1],'CH2_peak_divisions':divs[2],
                         'CH1_code_min':row['CH1_code_min'],'CH1_code_max':row['CH1_code_max'],
                         'CH2_code_min':row['CH2_code_min'],'CH2_code_max':row['CH2_code_max'],
                         'range_iterations':attempt})
    return row,scales[1],scales[2]

def max_sine_code(summary,ch):
    return max(abs(int(summary[f'CH{ch}_code_min'])),abs(int(summary[f'CH{ch}_code_max'])))

def valid_result(r):
    flags=[]
    clipped1=r['CH1_code_min']<=-32768 or r['CH1_code_max']>=32767
    clipped2=r['CH2_code_min']<=-32768 or r['CH2_code_max']>=32767
    if clipped1:flags.append('CH1 ADC clip')
    if clipped2:flags.append('CH2 ADC clip')
    if not np.isfinite(r['f0_measured_Hz']) or not 450000<=r['f0_measured_Hz']<=550000:flags.append('f0 outside scan span')
    if not np.isfinite(r['Hmax']) or r['Hmax']<=0:flags.append('invalid Hmax')
    if not np.isfinite(r['Vin_resonance_Vpk']) or r['Vin_resonance_Vpk']<=0:flags.append('invalid Vin')
    if not np.isfinite(r['Vout_resonance_Vpk']) or r['Vout_resonance_Vpk']<=0:flags.append('invalid Vout')
    if not np.isfinite(r['f_lower_3dB_Hz']) or not np.isfinite(r['f_upper_3dB_Hz']) or r['f_upper_3dB_Hz']<=r['f_lower_3dB_Hz']:
        flags.append('invalid -3 dB crossings')
    return not flags,'; '.join(flags)

def save_figures(df):
    order={1:'o',2:'s',3:'^'};colors={1:'C0',2:'C1',3:'C2'}
    grouped=df.groupby('amplitude_programmed_Vpp').agg(
        f0_mean=('f0_measured_Hz','mean'),f0_sd=('f0_measured_Hz','std'),
        h_mean=('Hmax','mean'),h_sd=('Hmax','std'),
        i_mean=('I_resonance_pk_A','mean'),i_sd=('I_resonance_pk_A','std'))
    fig,ax=plt.subplots(figsize=(8.5,5.3))
    for rep,g in df.groupby('repeat'):
        g=g.sort_values('amplitude_programmed_Vpp')
        ax.plot(g.amplitude_programmed_Vpp,g.f0_measured_Hz/1000,marker=order[int(rep)],ls='none',color=colors[int(rep)],label=f'Repeat {rep}')
    ax.errorbar(grouped.index,grouped.f0_mean/1000,yerr=grouped.f0_sd/1000,fmt='k-',marker='D',capsize=4,label='Mean ± sample SD')
    ax.set(xlabel='Programmed generator amplitude (Vpp)',ylabel='Phase-zero resonance f₀ (kHz)',title='Resonance frequency vs drive amplitude');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(OUT/'f0_vs_drive.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,5.3))
    for rep,g in df.groupby('repeat'):
        g=g.sort_values('amplitude_programmed_Vpp')
        ax.plot(g.amplitude_programmed_Vpp,g.Hmax,marker=order[int(rep)],ls='none',color=colors[int(rep)],label=f'Repeat {rep}')
    ax.errorbar(grouped.index,grouped.h_mean,yerr=grouped.h_sd,fmt='k-',marker='D',capsize=4,label='Mean ± sample SD')
    ax.set(xlabel='Programmed generator amplitude (Vpp)',ylabel='Hmax = |Vout/Vin|',title='Peak transfer vs drive amplitude');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(OUT/'Hmax_vs_drive.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,5.3))
    ax.errorbar(grouped.i_mean*1000,grouped.f0_mean/1000,xerr=grouped.i_sd*1000,yerr=grouped.f0_sd/1000,fmt='o-',capsize=4)
    for amp,(x,y) in zip(grouped.index,zip(grouped.i_mean*1000,grouped.f0_mean/1000)):
        ax.annotate(f'{amp:g} Vpp',(x,y),xytext=(4,5),textcoords='offset points',fontsize=8)
    ax.set(xlabel='Resonance current (mA peak; Vout,pk / 5.1 Ω)',ylabel='Phase-zero resonance f₀ (kHz)',title='Resonance frequency vs measured current');ax.grid(alpha=.25);fig.tight_layout();fig.savefig(OUT/'f0_vs_current.png',dpi=180);plt.close(fig)
    return grouped

def main():
    if OUT.exists():raise FileExistsError(f'Output directory already exists; refusing to overwrite: {OUT}')
    OUT.mkdir(parents=True)
    rm=scope=gen=None
    try:
        rm=pyvisa.ResourceManager();scope=rm.open_resource(method.base.SCOPE_ADDR);gen=rm.open_resource(method.base.GEN_ADDR)
        scope.timeout=120000;gen.timeout=15000;scope.chunk_size=4*1024*1024
        log('Fresh supporting measurement run; no existing result CSVs or measurement records loaded.')
        log('Scope '+scope.query('*IDN?').strip());log('Generator '+gen.query('*IDN?').strip())
        log('Initial generator: '+gen.query('C2:BSWV?').strip()+'; '+gen.query('C2:MDWV?').strip()+'; '+gen.query('C2:OUTP?').strip())
        scope.write(':CHAN1:COUP DC');scope.write(':CHAN2:COUP DC')
        scope.write(':CHAN1:PROB 1');scope.write(':CHAN2:PROB 1')
        scope.write(':TIM:SCAL 2E-4');scope.write(f':WAV:POIN {N}');scope.write(':WAV:INT 1');scope.write(':WAV:STAR 0')
        scope.write(':TRIG:RUN')
        scope.timeout=120000
        cfg={'timebase_s_div':float(scope.query(':TIM:SCAL?')),'sample_rate_Hz':float(scope.query(':ACQ:SRAT?')),
             'points':int(float(scope.query(':ACQ:POIN?'))),'wave_points':int(float(scope.query(':WAV:POIN?'))),
             'wave_interval':int(float(scope.query(':WAV:INT?'))),'CH1_impedance':scope.query(':CHAN1:IMP?').strip(),
             'CH2_impedance':scope.query(':CHAN2:IMP?').strip()}
        log('Scope acquisition setup '+json.dumps(cfg,sort_keys=True))
        if not all(('1M' in z.upper() or 'ONEMEG' in z.upper()) for z in (cfg['CH1_impedance'],cfg['CH2_impedance'])):
            raise RuntimeError(f'Scope channel input impedance is not confirmed as 1 MOhm: {cfg}')
        if not np.isclose(cfg['sample_rate_Hz'],FS,rtol=0,atol=1) or cfg['wave_points']!=N or cfg['wave_interval']!=1:
            raise RuntimeError(f'Scope did not accept 50 MSa/s, {N} points, interval 1: {cfg}')
        method.auto_range_at=fresh_auto_range_at
        rows=[];harm_rows=[];raw_rows=[];last_center=NOMINAL_F0
        sequence=[(amp,rep) for amp in AMPS for rep in (1,2,3)]
        for order,(amp,rep) in enumerate(sequence,1):
            log(f'START {order}/18: amplitude={amp:g} Vpp, repeat={rep}; nominal starting center={last_center:.2f} Hz')
            time.sleep(.2)
            row,harm,raw=method.one_sweep(scope,gen,amp,f'repeat_{rep}',order,last_center,RAW)
            resonance=next(r for r in raw if r['label']=='resonance_measurement')
            row.update({f'CH{ch}_code_{edge}':resonance[f'CH{ch}_code_{edge}']
                        for ch in (1,2) for edge in ('min','max')})
            scan_clips=[r for r in raw if r['CH1_code_min']<=-32768 or r['CH1_code_max']>=32767 or
                        r['CH2_code_min']<=-32768 or r['CH2_code_max']>=32767]
            accepted,reason=valid_result(row)
            if scan_clips:
                accepted=False;reason=(reason+'; ' if reason else '')+f'ADC endpoint clipping in {len(scan_clips)} frequency-scan captures'
            range_info=RANGE_CHECKS[-1]
            row.update({'repeat':rep,'scope_sample_rate_Hz':cfg['sample_rate_Hz'],'scope_record_points':N,
                        'scope_record_duration_s':N/cfg['sample_rate_Hz'],'clipped_CH1':row['CH1_code_min']<=-32768 or row['CH1_code_max']>=32767,
                        'clipped_CH2':row['CH2_code_min']<=-32768 or row['CH2_code_max']>=32767,
                        'measurement_valid':accepted,'validity_flags':reason,
                        'range_check_CH1_peak_divisions':range_info['CH1_peak_divisions'],
                        'range_check_CH2_peak_divisions':range_info['CH2_peak_divisions'],
                        'range_check_iterations':range_info['range_iterations'],
                        'frequency_scan_clipped_capture_count':len(scan_clips),
                        'Vout_resonance_Vpp':2*row['Vout_resonance_Vpk'],
                        'I_resonance_pk_A_from_Vout_over_R':row['Vout_resonance_Vpk']/R_EXT,
                        'I_resonance_rms_A_from_Vout_over_R':row['Vout_resonance_Vpk']/(R_EXT*np.sqrt(2)),
                        'Hmax_fit_valid':bool(row['complex_fit_reliable']),
                        'timestamp_utc':datetime.now(timezone.utc).isoformat(timespec='milliseconds')})
            harm['repeat']=rep;harm['measurement_valid']=accepted
            for r in raw:r['repeat']=rep
            rows.append(row);harm_rows.append(harm);raw_rows.extend(raw)
            last_center=float(row['f0_programmed_Hz'])
            log(f'COMPLETE {amp:g} Vpp repeat {rep}: f0={row["f0_measured_Hz"]:.2f} Hz (programmed {row["f0_programmed_Hz"]:.2f}), Hmax={row["Hmax"]:.6f}, Vin={row["Vin_resonance_Vpp"]:.6g} Vpp, Vout={row["Vout_resonance_Vpk"]:.6g} Vpk, Ipk={row["I_resonance_pk_A_from_Vout_over_R"]*1000:.4g} mA, scales={row["CH1_scale_V_div"]:g}/{row["CH2_scale_V_div"]:g} V/div, clip={row["clipped_CH1"]}/{row["clipped_CH2"]}, complex-fit-reliable={row["Hmax_fit_valid"]}; validity={reason or "PASS"}')
            pd.DataFrame(rows).to_csv(CSV,index=False)
            pd.DataFrame(raw_rows).to_csv(RAW,index=False)
            pd.DataFrame(harm_rows).to_csv(HARM,index=False)
        df=pd.DataFrame(rows)
        if len(df)!=len(AMPS)*REPEATS:raise RuntimeError(f'Incomplete repeats: {len(df)}')
        df.to_csv(CSV,index=False);pd.DataFrame(raw_rows).to_csv(RAW,index=False);pd.DataFrame(harm_rows).to_csv(HARM,index=False)
        grouped=save_figures(df)
        clips=bool(df[['clipped_CH1','clipped_CH2']].to_numpy().any())
        invalid=df[~df.measurement_valid]
        trend=float(grouped.loc[5.0,'f0_mean']-grouped.loc[0.2,'f0_mean'])
        trend_pct=100*trend/grouped.loc[0.2,'f0_mean']
        low_high_h=float(grouped.loc[5.0,'h_mean']-grouped.loc[0.2,'h_mean'])
        md=['# Fresh drive-dependence supporting measurement','',f'- Run completed: {datetime.now(timezone.utc).isoformat(timespec="seconds")}.','- Only this run’s acquisitions were used; no prior measurement result files were loaded.','- Method: existing pure-sine workflow with AM OFF, sine-fit complex transfer H, coarse frequency sweep, local phase-zero refinement, Hmax from refined scan, and direct CH1/CH2 resonance capture.','- Scope ranges were freshly auto-ranged from the current run’s captured waveform codes at each amplitude; per-repeat ranges are in the CSV.','',f'- Acquisitions: **{len(df)} / 18** ({len(AMPS)} amplitudes × {REPEATS} repeats).',f'- Validity: **{len(df)-len(invalid)} passed**, {len(invalid)} flagged; ADC clipping detected: **{"yes" if clips else "no"}**.',f'- Phase-zero f0 shift, 0.2 to 5 Vpp: **{trend/1000:+.4f} kHz** ({trend_pct:+.4f}%).',f'- Mean Hmax change, 0.2 to 5 Vpp: **{low_high_h:+.6g}**.',f'- Vout peak current calculated as Vout,pk / 5.1 Ω; both peak and RMS values are recorded.','', '| Drive (Vpp) | f0 mean (kHz) | f0 sample SD (Hz) | Hmax mean | Hmax sample SD | Vin at f0 mean (Vpp) | Vout at f0 mean (Vpk) | Ires mean (mA pk) | Scope scales CH1/CH2 (V/div) | Valid repeats |','|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|']
        for amp,g in df.groupby('amplitude_programmed_Vpp',sort=True):
            scales=', '.join(f'{a:g}/{b:g}' for a,b in sorted(set(zip(g.CH1_scale_V_div,g.CH2_scale_V_div))))
            md.append(f'| {amp:g} | {g.f0_measured_Hz.mean()/1000:.5f} | {g.f0_measured_Hz.std(ddof=1):.3g} | {g.Hmax.mean():.6f} | {g.Hmax.std(ddof=1):.3g} | {g.Vin_resonance_Vpp.mean():.6g} | {g.Vout_resonance_Vpk.mean():.6g} | {g.I_resonance_pk_A_from_Vout_over_R.mean()*1000:.5g} | {scales} | {int(g.measurement_valid.sum())}/3 |')
        md += ['', '## Quality checks','']
        if invalid.empty:md.append('- All repeated measurements passed finite-result, in-scan-range, valid Hmax and -3 dB crossing checks.')
        else:
            for _,r in invalid.iterrows():md.append(f'- Flagged {r.amplitude_programmed_Vpp:g} Vpp repeat {int(r.repeat)}: {r.validity_flags}.')
        if not clips:md.append('- No ADC endpoint clipping was found in the recorded resonance waveforms; code extrema and scope ranges are in the CSV.')
        md += ['', '## Files','', '- `drive_dependence_support.csv` — three fresh resonance results at each drive amplitude.','- `drive_dependence_frequency_scans.csv` — full per-frequency raw analysis rows from the phase, Hmax, and bandwidth scans.','- `drive_dependence_harmonics.csv` — resonance harmonic checks captured by the reused workflow.','- `f0_vs_drive.png`, `Hmax_vs_drive.png`, `f0_vs_current.png` — requested plots.','- `drive_dependence_support.log` — readbacks, scope settings, range optimization, and acquisition status.','']
        SUMMARY.write_text('\n'.join(md),encoding='utf-8')
        log(f'OUTPUTS complete: {CSV}, {RAW}, {SUMMARY}')
    finally:
        # Leave the apparatus at the last fresh run setting; do not restore an earlier state.
        try:
            if scope is not None:scope.write(':TRIG:RUN')
        except Exception:pass
        for inst in (gen,scope):
            if inst is not None:
                try:inst.close()
                except Exception:pass
        if rm is not None:
            try:rm.close()
            except Exception:pass
        if RUNLOG:LOG.write_text('\n'.join(RUNLOG)+'\n',encoding='utf-8')

if __name__=='__main__':main()
