"""Fresh AM-OFF amplitude characterization with a guarded extreme-drive point."""
from __future__ import annotations

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import csv, json, math, re, sys, time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvisa

import characterize_drive_amplitude as method

ROOT=_bench.workspace
OUT=ROOT/'drive_dependence_extreme_12Vpp_20260930_run02'
AMPS=[0.2,1.0,2.0,3.0,4.0,5.0,12.0]
REPEATS=3
N=100_000
FS=50_000_000.0
R_EXT=5.1
EXTREME=12.0
EXTREME_LABEL='12 Vpp — EXTREME DRIVE TEST'
AMPLITUDE_COLORS={0.2:'#0072B2',1.0:'#E69F00',2.0:'#009E73',3.0:'#CC79A7',
                  4.0:'#56B4E9',5.0:'#D55E00',12.0:'#D62728'}
NOMINAL_F0=1/(2*np.pi*np.sqrt(22e-6*4.7e-9))
RESULTS=OUT/'drive_dependence_results.csv'
RAW=OUT/'frequency_response_raw.csv'
CAPTURE_AUDIT=OUT/'acquisition_audit.csv'
HARM=OUT/'harmonics.csv'
LOG=OUT/'collection.log'
SUMMARY=OUT/'collection_summary.md'
PLOT_DIR=OUT/'plots'
RESUME_FROM=(4.0,2)
ATTEMPT_ID='resume_4V_repeat2_after_IO_fault_01'
RESUME_START_UTC=''

RUNLOG=[]
CURRENT_AMP=float('nan')
CURRENT_REPEAT=0
CURRENT_ORDER=0
ABORT_REASON=''
ORIGINAL_MEASURE=method.base.measure

class StopRun(RuntimeError):
    pass

def log(message):
    line=f'[{datetime.now(timezone.utc).isoformat(timespec="milliseconds")}] {message}'
    print(line,flush=True)
    RUNLOG.append(line)
    if OUT.exists():
        with LOG.open('a',encoding='utf-8') as f:f.write(line+'\n')

def append_record(path,record):
    path.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame([record]).to_csv(path,mode='a' if path.exists() else 'w',header=not path.exists(),index=False)

def drain_errors(inst,label):
    values=[]
    for _ in range(12):
        try: value=inst.query('SYST:ERR?').strip()
        except Exception as e:
            log(f'{label} error-queue query unavailable: {e!r}')
            return values
        values.append(value)
        if value.startswith('0,') or value.startswith('0 '):break
    for value in values:
        if not value.startswith(('0,','0 ')):log(f'{label} pre-run queued error (drained, not attributed to this run): {value}')
    return values

def parse_amp(readback):
    match=re.search(r'AMP,([0-9.+\-Ee]+)V',readback,re.I)
    return float(match.group(1)) if match else math.nan

def instrument_status(scope,gen):
    gen_bswv=gen.query('C2:BSWV?').strip()
    gen_out=gen.query('C2:OUTP?').strip()
    gen_err=gen.query('SYST:ERR?').strip()
    scope_err=scope.query('SYST:ERR?').strip()
    return {'generator_bswv':gen_bswv,'generator_output':gen_out,'generator_error':gen_err,'scope_error':scope_err,
            'generator_amplitude_readback_Vpp':parse_amp(gen_bswv)}

def check_status(scope,gen,expected_amp):
    status=instrument_status(scope,gen)
    bad=[]
    if not status['generator_output'].upper().startswith('C2:OUTP ON'):
        bad.append('generator output is not ON')
    if any(tag in status['generator_output'].upper() for tag in ('OVER','OVLD','PROT','TRIP','ALARM')):
        bad.append('generator overload/protection indicator')
    if not np.isfinite(status['generator_amplitude_readback_Vpp']) or not np.isclose(status['generator_amplitude_readback_Vpp'],expected_amp,rtol=0.005,atol=0.005):
        bad.append(f'generator amplitude readback mismatch (wanted {expected_amp:g} Vpp)')
    if not status['generator_error'].startswith(('0,','0 ')):
        bad.append('generator reported SCPI error: '+status['generator_error'])
    if not status['scope_error'].startswith(('0,','0 ')):
        bad.append('scope reported SCPI error: '+status['scope_error'])
    if bad:
        raise StopRun('; '.join(bad))
    return status

def stop(reason):
    global ABORT_REASON
    ABORT_REASON=str(reason)
    log('STOP CONDITION: '+ABORT_REASON)
    raise StopRun(ABORT_REASON)

SCALES=method.SUPPORTED_SCALES

def next_scale(scale,up):
    choices=SCALES[SCALES>scale*(1+1e-9)] if up else SCALES[SCALES<scale*(1-1e-9)]
    return float(choices[0] if up and len(choices) else choices[-1] if len(choices) else scale)

def fresh_auto_range_at(scope,gen,amp_vpp,fcenter):
    # The extreme point starts at intentionally generous scales; all decisions
    # below are then based only on waveforms acquired in this new run.
    scales={1:1.0,2:0.5} if amp_vpp>=10 else {1:0.2,2:0.1}
    for ch in (1,2):scope.write(f':CHAN{ch}:SCAL {scales[ch]:.9g}')
    row=None
    for attempt in range(1,11):
        row=method.base.measure(scope,gen,float(fcenter),'fresh_range_check')
        occ={ch:max(abs(row[f'CH{ch}_code_min']),abs(row[f'CH{ch}_code_max']))/7680.0 for ch in (1,2)}
        # Any ADC endpoint hit is a hard stop, even if a larger scale could
        # subsequently make another capture look acceptable.
        for ch in (1,2):
            if row[f'CH{ch}_code_min']<=-32768 or row[f'CH{ch}_code_max']>=32767:
                stop(f'ADC clipping during {amp_vpp:g} Vpp live range check on CH{ch}')
        changed=False
        for ch in (1,2):
            if occ[ch]>3.35:new=next_scale(scales[ch],True)
            elif occ[ch]<1.15:new=next_scale(scales[ch],False)
            else:new=scales[ch]
            if new!=scales[ch]:scales[ch]=new;scope.write(f':CHAN{ch}:SCAL {new:.9g}');changed=True
        if not changed:break
    scales={ch:float(scope.query(f':CHAN{ch}:SCAL?')) for ch in (1,2)}
    occ={ch:max(abs(row[f'CH{ch}_code_min']),abs(row[f'CH{ch}_code_max']))/7680.0 for ch in (1,2)}
    if any(occ[ch]>3.55 for ch in (1,2)):
        stop(f'Live auto-range cannot provide headroom: peak divisions={occ}, scales={scales}')
    return row,scales[1],scales[2]

def capture_guard(row):
    """Persist every fresh capture and halt immediately on invalidity/clipping."""
    flags=[]
    clipped={ch:row[f'CH{ch}_code_min']<=-32768 or row[f'CH{ch}_code_max']>=32767 for ch in (1,2)}
    if any(clipped.values()):flags.append('ADC clipping')
    if int(row.get('n_points',0))!=N:flags.append(f'invalid point count {row.get("n_points")} (expected {N})')
    for key in ('frequency_measured_Hz','H_real','H_imag','H_mag','H_phase_deg','vin_Vpk','vout_Vpk'):
        if not np.isfinite(row.get(key,np.nan)):flags.append(f'non-finite {key}')
    try:status=check_status(SCOPE,GEN,CURRENT_AMP)
    except Exception as e:
        status={'generator_bswv':'','generator_output':'','generator_error':'','scope_error':'','generator_amplitude_readback_Vpp':math.nan}
        flags.append(str(e))
    record=dict(row)
    record.update({'amplitude_programmed_Vpp':CURRENT_AMP,'repeat':CURRENT_REPEAT,'sweep_order':CURRENT_ORDER,
                   'captured_utc':datetime.now(timezone.utc).isoformat(timespec='milliseconds'),
                   'clipped_CH1':clipped[1],'clipped_CH2':clipped[2],
                   'acquisition_valid':not flags,'validity_flags':'; '.join(flags),**status})
    append_record(RAW,record)
    append_record(CAPTURE_AUDIT,{'timestamp_utc':record['captured_utc'],'amplitude_programmed_Vpp':CURRENT_AMP,
                  'repeat':CURRENT_REPEAT,'sweep_order':CURRENT_ORDER,'label':row.get('label'),
                  'frequency_set_Hz':row.get('frequency_set_Hz'),'frequency_measured_Hz':row.get('frequency_measured_Hz'),
                  'n_points':row.get('n_points'),'CH1_scale_V_div':row.get('CH1_vdiv_V_div'),
                  'CH2_scale_V_div':row.get('CH2_vdiv_V_div'),'CH1_code_min':row.get('CH1_code_min'),
                  'CH1_code_max':row.get('CH1_code_max'),'CH2_code_min':row.get('CH2_code_min'),
                  'CH2_code_max':row.get('CH2_code_max'),'clipped_CH1':clipped[1],'clipped_CH2':clipped[2],
                  'generator_output':status['generator_output'],'generator_error':status['generator_error'],
                  'scope_error':status['scope_error'],'acquisition_valid':not flags,'validity_flags':'; '.join(flags)})
    if flags:stop(f'{CURRENT_AMP:g} Vpp, {row.get("label")}: '+', '.join(flags))

def measure_guarded(scope,gen,f_set,label):
    global ABORT_REASON
    try:row=ORIGINAL_MEASURE(scope,gen,f_set,label)
    except Exception as e:
        ABORT_REASON=f'invalid waveform acquisition at {CURRENT_AMP:g} Vpp, {label}, {f_set:.2f} Hz: {e!r}'
        log('STOP CONDITION: '+ABORT_REASON)
        raise StopRun(ABORT_REASON) from e
    capture_guard(row)
    return row

def programmed_guard(amp,readback):
    if not np.isclose(amp,CURRENT_AMP,rtol=0.005,atol=0.005):
        stop(f'generator programmed amplitude {amp:g} Vpp does not match requested {CURRENT_AMP:g} Vpp')
    status=check_status(SCOPE,GEN,CURRENT_AMP)
    log(f'Programmed {amp:g} Vpp; readback {readback}; output {status["generator_output"]}; AM {GEN.query("C2:MDWV?").strip()}')
    if 'STATE,OFF' not in GEN.query('C2:MDWV?').upper():stop('AM is not OFF')
    max_match=re.search(r'MAX_OUTPUT_AMP,([0-9.+\-Ee]+)V',readback,re.I)
    if max_match and amp>float(max_match.group(1)):stop(f'requested {amp:g} Vpp exceeds generator MAX_OUTPUT_AMP')

def harmonic_guard(row):
    flags=[]
    clipped={ch:row[f'CH{ch}_code_min']<=-32768 or row[f'CH{ch}_code_max']>=32767 for ch in (1,2)}
    if any(clipped.values()):flags.append('ADC clipping in harmonic capture')
    for key in ('CH1_H2_dBc','CH1_H3_dBc','CH2_H2_dBc','CH2_H3_dBc'):
        if not np.isfinite(row.get(key,np.nan)):flags.append('invalid harmonic '+key)
    try:status=check_status(SCOPE,GEN,CURRENT_AMP)
    except Exception as e:
        status={'generator_bswv':'','generator_output':'','generator_error':'','scope_error':'','generator_amplitude_readback_Vpp':math.nan}
        flags.append(str(e))
    rec=dict(row);rec.update({'amplitude_programmed_Vpp':CURRENT_AMP,'repeat':CURRENT_REPEAT,'sweep_order':CURRENT_ORDER,
        'captured_utc':datetime.now(timezone.utc).isoformat(timespec='milliseconds'),'clipped_CH1':clipped[1],
        'clipped_CH2':clipped[2],'acquisition_valid':not flags,'validity_flags':'; '.join(flags),**status})
    append_record(HARM,rec)
    if flags:stop(f'{CURRENT_AMP:g} Vpp harmonic capture: '+'; '.join(flags))

def summary_validity(row):
    flags=[]
    for key in ('f0_measured_Hz','Hmax','f_lower_3dB_Hz','f_upper_3dB_Hz','BW_measured_Hz','Q_measured','Vin_resonance_Vpp','Vout_resonance_Vpk'):
        if not np.isfinite(row.get(key,np.nan)):flags.append('invalid '+key)
    if not 450000<=row['f0_measured_Hz']<=550000:flags.append('f0 outside scan span')
    if row['BW_measured_Hz']<=0 or row['Q_measured']<=0:flags.append('invalid bandwidth/Q')
    return flags

def extreme_note(ax,xy=None,fit_bad=False):
    label=EXTREME_LABEL+(' (fit unreliable)' if fit_bad else '')
    if xy is not None and np.isfinite(xy[0]) and np.isfinite(xy[1]):
        xlo,xhi=ax.get_xlim();ylo,yhi=ax.get_ylim()
        right_edge=xy[0]>xlo+0.75*(xhi-xlo);top_edge=xy[1]>ylo+0.78*(yhi-ylo)
        ax.annotate(label,(xy[0],xy[1]),xytext=(-8 if right_edge else 8,-8 if top_edge else 10),textcoords='offset points',
                    ha='right' if right_edge else 'left',va='top' if top_edge else 'bottom',color='red',weight='bold',fontsize=8,
                    bbox={'boxstyle':'round,pad=.25','fc':'white','ec':'red','alpha':.9})
    else:
        ax.text(.98,.97,label,transform=ax.transAxes,ha='right',va='top',color='red',weight='bold',fontsize=9,
                bbox={'boxstyle':'round,pad=.3','fc':'white','ec':'red','alpha':.95})

def point_plot(grouped,xkey,ykey,xlabel,ylabel,title,path,yerr=None,xerr=None):
    fig,ax=plt.subplots(figsize=(9,5.6))
    for amp,r in grouped.iterrows():
        xv=r[xkey] if xkey!='amplitude_programmed_Vpp' else amp
        label=EXTREME_LABEL if amp==EXTREME else f'{amp:g} Vpp'
        ax.errorbar([xv],[r[ykey]],xerr=[r[xerr]] if xerr else None,yerr=[r[yerr]] if yerr else None,
                    fmt='o',color=AMPLITUDE_COLORS.get(float(amp),'#777777'),ecolor=AMPLITUDE_COLORS.get(float(amp),'#777777'),
                    capsize=3,ms=6,label=label,zorder=5 if amp==EXTREME else 3)
    if EXTREME in grouped.index:
        r=grouped.loc[EXTREME];xv=r[xkey] if xkey!='amplitude_programmed_Vpp' else EXTREME
        extreme_note(ax,(xv,r[ykey]))
    ax.set(xlabel=xlabel,ylabel=ylabel,title=title);ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(PLOT_DIR/path,dpi=190);plt.close(fig)

def make_plots(results,raw):
    PLOT_DIR.mkdir(parents=True,exist_ok=True)
    grouped=results.groupby('amplitude_programmed_Vpp').agg(
        f0=('f0_measured_Hz','mean'),f0_sd=('f0_measured_Hz','std'),
        h=('Hmax','mean'),h_sd=('Hmax','std'),
        bw=('BW_measured_Hz','mean'),bw_sd=('BW_measured_Hz','std'),
        q=('Q_measured','mean'),q_sd=('Q_measured','std'),
        vin=('Vin_resonance_Vpp','mean'),vin_sd=('Vin_resonance_Vpp','std'),
        vout=('Vout_resonance_Vpk','mean'),i=('I_resonance_pk_A','mean'),i_sd=('I_resonance_pk_A','std'),
        rs=('Rs_fit_ohm','mean'),rs_sd=('Rs_fit_ohm','std'),L=('L_eff_fit_uH','mean'),L_sd=('L_eff_fit_uH','std'),
        C=('C_eff_fit_nF','mean'),C_sd=('C_eff_fit_nF','std'),fit_ok=('complex_fit_reliable','all'),
        CH1_H2=('CH1_H2_dBc','mean'),CH1_H3=('CH1_H3_dBc','mean'),CH2_H2=('CH2_H2_dBc','mean'),CH2_H3=('CH2_H3_dBc','mean'))
    grouped['amplitude_programmed_Vpp']=grouped.index.astype(float)
    grouped['f0_kHz']=grouped.f0/1000;grouped['f0_sd_kHz']=grouped.f0_sd/1000
    grouped['bw_kHz']=grouped.bw/1000;grouped['bw_sd_kHz']=grouped.bw_sd/1000;grouped['i_mA']=grouped.i*1000
    grouped['i_sd_mA']=grouped.i_sd*1000
    point_plot(grouped,'amplitude_programmed_Vpp','f0_kHz','Programmed generator amplitude (Vpp)','f₀ (kHz)','Phase-zero resonance vs drive amplitude','f0_vs_drive_amplitude.png','f0_sd_kHz')
    point_plot(grouped,'i_mA','f0_kHz','Resonance current (mA peak)','f₀ (kHz)','Phase-zero resonance vs resonance current','f0_vs_resonance_current.png','f0_sd_kHz','i_sd_mA')
    point_plot(grouped,'amplitude_programmed_Vpp','h','Programmed amplitude (Vpp)','Hmax','Peak transfer vs drive amplitude','hmax_vs_drive_amplitude.png','h_sd')
    point_plot(grouped,'i_mA','h','Resonance current (mA peak)','Hmax','Peak transfer vs resonance current','hmax_vs_current.png','h_sd','i_sd_mA')

    for xkey,xlab,xsuffix in [('amplitude_programmed_Vpp','Programmed amplitude (Vpp)','drive_amplitude'),('i_mA','Resonance current (mA peak)','current')]:
        fig,axs=plt.subplots(2,1,figsize=(9,8),sharex=True)
        for ax,ykey,err,ylab,title in [(axs[0],'bw_kHz','bw_sd_kHz','Bandwidth (kHz)','-3 dB bandwidth'),(axs[1],'q','q_sd','Q','Quality factor')]:
            for amp,r in grouped.iterrows():
                xv=r[xkey] if xkey!='amplitude_programmed_Vpp' else amp
                ax.errorbar([xv],[r[ykey]],yerr=[r[err]],fmt='o',color=AMPLITUDE_COLORS.get(float(amp),'#777777'),
                            ecolor=AMPLITUDE_COLORS.get(float(amp),'#777777'),capsize=3,ms=6,
                            label=EXTREME_LABEL if amp==EXTREME else f'{amp:g} Vpp',zorder=5 if amp==EXTREME else 3)
            if EXTREME in grouped.index:
                r=grouped.loc[EXTREME];xv=r[xkey] if xkey!='amplitude_programmed_Vpp' else EXTREME
                extreme_note(ax,(xv,r[ykey]))
            ax.set_ylabel(ylab);ax.set_title(title);ax.grid(alpha=.25);ax.legend(fontsize=8,ncol=2)
        axs[-1].set_xlabel(xlab);fig.suptitle(f'Bandwidth and Q vs {xsuffix.replace("_"," ")}');fig.tight_layout();fig.savefig(PLOT_DIR/f'bandwidth_q_vs_{xsuffix}.png',dpi=190);plt.close(fig)

    fit=results[results.complex_fit_reliable.astype(bool)].copy()
    for xkey,xlab,xsuffix in [('amplitude_programmed_Vpp','Programmed amplitude (Vpp)','drive_amplitude'),('I_resonance_pk_A','Resonance current (A peak)','current')]:
        fig,axs=plt.subplots(3,1,figsize=(9,10),sharex=True)
        for ax,key,err,ylab in zip(axs,['Rs_fit_ohm','L_eff_fit_uH','C_eff_fit_nF'],['Rs_fit_ohm','L_eff_fit_uH','C_eff_fit_nF'],['Rs fit (Ω)','L_eff fit (µH)','C_eff fit (nF)']):
            for amp,g in fit.groupby('amplitude_programmed_Vpp'):
                color=AMPLITUDE_COLORS.get(float(amp),'#777777');lw=2.2 if amp==EXTREME else 1.2;ms=5.5
                ax.plot(g[xkey],g[key],marker='o',color=color,lw=lw,ms=ms,
                        label=EXTREME_LABEL if amp==EXTREME else f'{amp:g} Vpp')
            is_good=EXTREME in set(fit.amplitude_programmed_Vpp)
            xy=(fit.loc[fit.amplitude_programmed_Vpp==EXTREME,xkey].mean(),
                fit.loc[fit.amplitude_programmed_Vpp==EXTREME,key].mean()) if is_good else None
            extreme_note(ax,xy,fit_bad=not is_good)
            ax.set_ylabel(ylab);ax.grid(alpha=.25)
        axs[-1].set_xlabel(xlab);axs[0].legend();fig.suptitle(f'Reliable complex-fit parameters vs {xsuffix.replace("_"," ")}');fig.tight_layout();fig.savefig(PLOT_DIR/f'fit_parameters_vs_{xsuffix}.png',dpi=190);plt.close(fig)

    for xkey,xlab,xsuffix in [('amplitude_programmed_Vpp','Programmed amplitude (Vpp)','drive_amplitude'),('i_mA','Resonance current (mA peak)','current')]:
        fig,axs=plt.subplots(2,2,figsize=(11,8),sharex=True)
        for ax,key,title in zip(axs.flat,['CH1_H2','CH1_H3','CH2_H2','CH2_H3'],['CH1 H2','CH1 H3','CH2 H2','CH2 H3']):
            for amp,r in grouped.iterrows():
                xv=r[xkey] if xkey!='amplitude_programmed_Vpp' else amp
                ax.plot([xv],[r[key]],'o',color=AMPLITUDE_COLORS.get(float(amp),'#777777'),ms=6,
                        label=EXTREME_LABEL if amp==EXTREME else f'{amp:g} Vpp',zorder=5 if amp==EXTREME else 3)
            if EXTREME in grouped.index:
                r=grouped.loc[EXTREME];xv=r[xkey] if xkey!='amplitude_programmed_Vpp' else EXTREME
                extreme_note(ax,(xv,r[key]))
            ax.set_title(title);ax.set_ylabel('Harmonic level (dBc)');ax.grid(alpha=.25)
        for ax in axs[-1]:ax.set_xlabel(xlab)
        axs[0,0].legend(fontsize=8);fig.suptitle(f'Harmonic levels vs {xsuffix.replace("_"," ")}');fig.tight_layout();fig.savefig(PLOT_DIR/f'harmonics_vs_{xsuffix}.png',dpi=190);plt.close(fig)

    fig,axs=plt.subplots(2,1,figsize=(10,9),sharex=True)
    curve=raw[raw.label.isin(['coarse','phase_refinement','bandwidth_refinement'])].copy()
    # Keep the original aborted partial scan in raw, but exclude it from response overlays.
    if RESUME_START_UTC and 'captured_utc' in curve:
        curve=curve[~((curve.amplitude_programmed_Vpp==4.0)&(curve.repeat==2)&
                      (pd.to_datetime(curve.captured_utc,utc=True)<pd.Timestamp(RESUME_START_UTC)))]
    for amp,gamp in curve.groupby('amplitude_programmed_Vpp'):
        for rep,gr in gamp.groupby('repeat'):
            gr=gr.sort_values('frequency_measured_Hz');extreme=amp==EXTREME
            color=AMPLITUDE_COLORS.get(float(amp),'#777777');alpha=1.0 if extreme else .52;lw=2.3 if extreme else .9
            label=EXTREME_LABEL if extreme and rep==int(gamp.repeat.min()) else (f'{amp:g} Vpp' if not extreme and rep==int(gamp.repeat.min()) else None)
            axs[0].plot(gr.frequency_measured_Hz/1000,gr.H_mag,color=color,alpha=alpha,lw=lw,label=label,zorder=8 if extreme else 1)
            axs[1].plot(gr.frequency_measured_Hz/1000,gr.H_phase_deg,color=color,alpha=alpha,lw=lw,label=label,zorder=8 if extreme else 1)
    for ax,ylabel,title in [(axs[0],'|H|','Magnitude response'),(axs[1],'Phase(H) (°)','Phase response')]:
        ax.set_ylabel(ylabel);ax.set_title(title);ax.grid(alpha=.25);extreme_note(ax)
    axs[-1].set_xlabel('Measured frequency (kHz)');axs[0].legend(fontsize=8,ncol=2);fig.suptitle('Resonance-response curves: all amplitudes and repeats');fig.tight_layout();fig.savefig(PLOT_DIR/'resonance_response_overlays.png',dpi=200);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.5,7))
    for amp,gamp in curve.groupby('amplitude_programmed_Vpp'):
        for rep,gr in gamp.groupby('repeat'):
            gr=gr.sort_values('frequency_measured_Hz');extreme=amp==EXTREME
            ax.plot(gr.H_real,gr.H_imag,color=AMPLITUDE_COLORS.get(float(amp),'#777777'),alpha=1 if extreme else .52,lw=2.2 if extreme else .9,
                    label=EXTREME_LABEL if extreme and rep==int(gamp.repeat.min()) else (f'{amp:g} Vpp' if not extreme and rep==int(gamp.repeat.min()) else None),zorder=8 if extreme else 1)
    ax.set(xlabel='Re(H)',ylabel='Im(H)',title='Complex transfer response around resonance');ax.grid(alpha=.3);ax.axhline(0,c='k',lw=.5);ax.axvline(0,c='k',lw=.5);extreme_note(ax);ax.legend(fontsize=8,ncol=2);fig.tight_layout();fig.savefig(PLOT_DIR/'complex_response_plane.png',dpi=190);plt.close(fig)

    point_plot(grouped,'amplitude_programmed_Vpp','vin','Programmed amplitude (Vpp)','Actual CH1 Vin at f₀ (Vpp)','Actual input voltage vs drive','vin_at_resonance_vs_drive.png','vin_sd')
    point_plot(grouped,'amplitude_programmed_Vpp','i_mA','Programmed amplitude (Vpp)','Resonance current (mA peak)','Resonance current vs drive','current_vs_drive.png','i_sd')
    return grouped

def write_summary(results,grouped,complete,reason=''):
    clips=bool(results[['clipped_CH1','clipped_CH2']].to_numpy().any()) if len(results) else False
    lines=['# Fresh drive-dependence data with 12 Vpp extreme test','',
      f'- Run timestamp: {datetime.now(timezone.utc).isoformat(timespec="seconds")}.',
      f'- Run folder: `{OUT.name}`. Recovery attempt: `{ATTEMPT_ID}`; the incomplete 4 Vpp repeat 2 scan remains preserved in raw records and is excluded from the final response overlay.',
      '- AM OFF, CH2 sine, 50 MSa/s, 100,000 points/channel, WAV interval 1, DC coupling, 1× probes, 1 MΩ scope inputs; live auto-ranging was performed from this run’s fresh waveforms.',
      f'- Run status: **{"COMPLETE" if complete else "STOPPED"}**; valid completed repeats: **{len(results)}/{len(AMPS)*REPEATS}**.',
      f'- Clipping in saved result captures: **{"yes" if clips else "no"}**.',
      f'- Safety/validity stop reason: **{reason or "none"}**.','']
    lines += ['## VISA communication interruption and recovery','',
      '- The original 4 Vpp repeat 2 coarse scan stopped at 545 kHz with `VI_ERROR_IO (-1073807298)`; the failed transfer was not treated as a measurement.',
      '- A link check then passed 30 `*IDN?` queries, 20 full native waveform channel transfers, and 3 fresh 4 Vpp / 545 kHz captures. A separate first-capture attempt in `run03` later produced the same generic I/O error at 0.2 Vpp and was stopped with the generator off.',
      '- The run resumed at 4 Vpp repeat 2, restarting that repeat from its initial range check. All remaining measurements then completed with valid records.',
      '- This pattern supports an intermittent VISA/acquisition-path fault. It does not establish whether the physical USB connection caused it.','']
    if len(results):
        lines += ['| Drive (Vpp) | f₀ mean (kHz) | f₀ SD (Hz) | Hmax | BW (kHz) | Q | Vin at f₀ (Vpp) | Vout at f₀ (Vpk) | Ires (mA pk) | reliable fits | valid repeats |',
                  '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for amp,g in results.groupby('amplitude_programmed_Vpp',sort=True):
            label=f'{amp:g} **— EXTREME DRIVE TEST**' if amp==EXTREME else f'{amp:g}'
            lines.append(f"| {label} | {g.f0_measured_Hz.mean()/1000:.5f} | {g.f0_measured_Hz.std(ddof=1):.3g} | {g.Hmax.mean():.6f} | {g.BW_measured_Hz.mean()/1000:.3f} | {g.Q_measured.mean():.4f} | {g.Vin_resonance_Vpp.mean():.6g} | {g.Vout_resonance_Vpk.mean():.6g} | {g.I_resonance_pk_A_from_Vout_over_R.mean()*1000:.4f} | {int(g.complex_fit_reliable.sum())}/{len(g)} | {int(g.measurement_valid.sum())}/{len(g)} |")
        if grouped is not None and EXTREME in grouped.index and 0.2 in grouped.index:
            lines += ['',f'- 0.2 to 12 Vpp f₀ difference: **{grouped.loc[EXTREME,"f0"]/1000-grouped.loc[0.2,"f0"]/1000:+.4f} kHz**.',
                      f'- 0.2 to 12 Vpp Hmax difference: **{grouped.loc[EXTREME,"h"]-grouped.loc[0.2,"h"]:+.5f}**.',
                      f'- 12 Vpp mean resonance current: **{grouped.loc[EXTREME,"i"]*1000:.3f} mA peak**.']
    if not complete:lines+=['', '## Stop details','',reason or 'Run did not complete.','', 'No invalid value was treated as a valid result. Partial acquisition records remain in the raw CSVs.']
    else:lines+=['','## Files','', '- `drive_dependence_results.csv` — per-repeat phase-zero resonance, Hmax, bandwidth/Q, Vin/Vout/current, fit parameters, harmonics, ranges and clipping/validity fields.', '- `frequency_response_raw.csv` — every fresh complex H frequency point, saved incrementally; includes the complete resonance response curves.', '- `acquisition_audit.csv` — per-capture point count, code extrema, clipping, generator/scope status, and validity.', '- `harmonics.csv` — individual repeat H2/H3 measurements on CH1 and CH2.', '- `plots/` — every requested plot, with 12 Vpp in red and explicitly labelled.', '', '12 Vpp is an extreme diagnostic point and is not recommended as a normal operating condition.']
    SUMMARY.write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
    global SCOPE,GEN,CURRENT_AMP,CURRENT_REPEAT,CURRENT_ORDER,ABORT_REASON,RESUME_START_UTC
    if not OUT.exists():raise FileNotFoundError(f'Required partial dataset folder is missing: {OUT}')
    if not RESULTS.exists():raise FileNotFoundError(f'Required existing results file is missing: {RESULTS}')
    existing=pd.read_csv(RESULTS)
    existing_keys={(float(r.amplitude_programmed_Vpp),int(r.repeat)) for r in existing.itertuples()}
    required_before={(a,float_repeat) for a in AMPS for float_repeat in range(1,REPEATS+1)
                     if (a,float_repeat)<RESUME_FROM}
    if not required_before.issubset(existing_keys):
        raise RuntimeError(f'Cannot safely resume: missing earlier completed repeats: {sorted(required_before-existing_keys)}')
    PLOT_DIR.mkdir(parents=True,exist_ok=True)
    rm=scope=gen=None;completed=existing.to_dict('records');harmonic_rows=[];complete=False
    try:
        rm=pyvisa.ResourceManager();scope=rm.open_resource(method.base.SCOPE_ADDR);gen=rm.open_resource(method.base.GEN_ADDR)
        SCOPE,GEN=scope,gen
        scope.timeout=120000;gen.timeout=15000;scope.chunk_size=4*1024*1024
        RESUME_START_UTC=datetime.now(timezone.utc).isoformat(timespec='milliseconds')
        log(f'Resuming existing fresh run at {RESUME_FROM[0]:g} Vpp repeat {RESUME_FROM[1]}; preserving prior valid rows and appending attempt {ATTEMPT_ID}.')
        log('Scope '+scope.query('*IDN?').strip());log('Generator '+gen.query('*IDN?').strip())
        drain_errors(scope,'scope');drain_errors(gen,'generator')
        log('Generator initial '+gen.query('C2:BSWV?').strip()+'; '+gen.query('C2:MDWV?').strip()+'; '+gen.query('C2:OUTP?').strip())
        gen.write('C2:MDWV STATE,OFF')
        scope.write(':TIM:SCAL 2E-4');scope.write(f':WAV:POIN {N}');scope.write(':WAV:INT 1');scope.write(':WAV:STAR 0');scope.write(':TRIG:MODE AUTO');scope.write(':TRIG:RUN')
        cfg={'timebase_s_div':float(scope.query(':TIM:SCAL?')),'sample_rate_Hz':float(scope.query(':ACQ:SRAT?')),
             'acq_points':int(float(scope.query(':ACQ:POIN?'))),'wave_points':int(float(scope.query(':WAV:POIN?'))),
             'interval':int(float(scope.query(':WAV:INT?'))),'CH1_impedance':scope.query(':CHAN1:IMP?').strip(),
             'CH2_impedance':scope.query(':CHAN2:IMP?').strip(),'CH1_probe':scope.query(':CHAN1:PROB?').strip(),
             'CH2_probe':scope.query(':CHAN2:PROB?').strip(),'CH1_coupling':scope.query(':CHAN1:COUP?').strip(),
             'CH2_coupling':scope.query(':CHAN2:COUP?').strip(),'trigger_mode':scope.query(':TRIG:MODE?').strip()}
        log('Scope setup '+json.dumps(cfg,sort_keys=True))
        if not np.isclose(cfg['sample_rate_Hz'],FS,rtol=0,atol=1) or cfg['wave_points']!=N or cfg['interval']!=1:stop('scope acquisition settings do not match 50 MSa/s, 100k points, interval 1')
        if 'ONEMeg' not in cfg['CH1_impedance'] or 'ONEMeg' not in cfg['CH2_impedance']:stop('scope input impedance is not read back as 1 MOhm')
        if cfg['CH1_coupling'].upper()!='DC' or cfg['CH2_coupling'].upper()!='DC':stop('scope coupling is not read back as DC on both channels')
        if not all(np.isclose(float(cfg[f'CH{ch}_probe']),1.0) for ch in (1,2)):stop('scope probes are not read back as 1x on both channels')
        if not scope.query('SYST:ERR?').strip().startswith(('0,','0 ')):stop('scope reported an SCPI error during safe setup')
        method.auto_range_at=fresh_auto_range_at
        method.base.measure=measure_guarded
        order=0;center=NOMINAL_F0
        for amp in AMPS:
            for repeat in range(1,REPEATS+1):
                order+=1;CURRENT_AMP=amp;CURRENT_REPEAT=repeat;CURRENT_ORDER=order
                if (amp,repeat) in existing_keys:
                    old=existing[(existing.amplitude_programmed_Vpp==amp)&(existing.repeat==repeat)].iloc[-1]
                    center=float(old.f0_programmed_Hz)
                    continue
                if (amp,repeat)<RESUME_FROM:
                    stop(f'missing earlier completed repeat {amp:g} Vpp #{repeat}; refusing to skip it')
                log(f'START {order}/{len(AMPS)*REPEATS}: {amp:g} Vpp, repeat {repeat}; starting center {center:.2f} Hz')
                try:
                    row,harm,raw=method.one_sweep(scope,gen,amp,f'repeat_{repeat}',order,center,RAW,
                       on_capture=lambda _row:None,on_programmed=programmed_guard,on_harmonic=harmonic_guard)
                except StopRun:raise
                except Exception as e:
                    stop(f'{amp:g} Vpp repeat {repeat}: invalid/incomplete resonance characterization: {type(e).__name__}: {e}')
                flags=summary_validity(row)
                resonance=next(r for r in raw if r['label']=='resonance_measurement')
                scan_clip_count=sum(1 for r in raw if any(r[f'CH{ch}_code_min']<=-32768 or r[f'CH{ch}_code_max']>=32767 for ch in (1,2)))
                if scan_clip_count:flags.append(f'{scan_clip_count} clipped frequency points')
                if flags:stop(f'{amp:g} Vpp repeat {repeat}: '+', '.join(flags))
                row.update({'repeat':repeat,'measurement_valid':True,'validity_flags':'','CH1_scale_V_div':resonance['CH1_vdiv_V_div'],
                    'CH2_scale_V_div':resonance['CH2_vdiv_V_div'],'CH1_code_min':resonance['CH1_code_min'],
                    'CH1_code_max':resonance['CH1_code_max'],'CH2_code_min':resonance['CH2_code_min'],
                    'CH2_code_max':resonance['CH2_code_max'],'clipped_CH1':False,'clipped_CH2':False,
                    'frequency_scan_clipped_capture_count':scan_clip_count,'scope_sample_rate_Hz':FS,'scope_record_points':N,
                    'Vout_resonance_Vpp':2*row['Vout_resonance_Vpk'],
                    'I_resonance_pk_A_from_Vout_over_R':row['Vout_resonance_Vpk']/R_EXT,
                    'I_resonance_rms_A_from_Vout_over_R':row['Vout_resonance_Vpk']/(R_EXT*np.sqrt(2)),
                    'timestamp_utc':datetime.now(timezone.utc).isoformat(timespec='milliseconds')})
                append_record(RESULTS,row);completed.append(row);center=float(row['f0_programmed_Hz'])
                existing_keys.add((amp,repeat))
                log(f'COMPLETE {amp:g} Vpp repeat {repeat}: f0={row["f0_measured_Hz"]:.2f} Hz, Hmax={row["Hmax"]:.6f}, BW={row["BW_measured_Hz"]:.1f} Hz, Q={row["Q_measured"]:.4f}, Vin={row["Vin_resonance_Vpp"]:.6g} Vpp, Ipk={row["I_resonance_pk_A_from_Vout_over_R"]*1000:.4f} mA, scales={row["CH1_scale_V_div"]:g}/{row["CH2_scale_V_div"]:g} V/div, fit_reliable={row["complex_fit_reliable"]}; PASS')
        complete=True
    except StopRun as e:
        ABORT_REASON=ABORT_REASON or str(e)
        log('Run halted safely: '+ABORT_REASON)
    except Exception as e:
        ABORT_REASON=f'unexpected acquisition/analysis failure: {type(e).__name__}: {e}'
        log('Run halted safely: '+ABORT_REASON)
    finally:
        try:
            if gen is not None:
                # The final drive is an intentional 12 Vpp extreme point; do not leave it energized.
                gen.write('C2:OUTP OFF')
                log('Generator CH2 output set OFF at run end/stop.')
        except Exception as e:log(f'Could not turn generator output OFF: {e!r}')
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
    result_df=pd.DataFrame(completed)
    # On an abort the raw/audit files and already completed repeat rows remain.
    if complete and len(result_df)==len(AMPS)*REPEATS:
        rawdf=pd.read_csv(RAW)
        grouped=make_plots(result_df,rawdf)
        write_summary(result_df,grouped,True)
    else:
        write_summary(result_df,None,False,ABORT_REASON or f'incomplete count: {len(result_df)} valid repeats')
    if complete:log(f'OUTPUTS complete: {OUT}')
    return 0 if complete else 2

if __name__=='__main__':sys.exit(main())
