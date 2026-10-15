"""Recompute formal AM statistics and drive summaries without hardware or waveforms."""
import argparse
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FORMAL = {'am_initial': (75,15), 'am_extended': (155,31)}
SIDEBANDS = ('T_lower','T_upper','T_mean')


def assert_close(actual, expected, name):
    if not np.allclose(actual, expected, rtol=1e-9, atol=1e-12, equal_nan=True):
        raise ValueError(f'Measured-table consistency check failed: {name}')


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def formal_summary(raw):
    rows=[]
    for fm,g in raw.groupby('fm_Hz',sort=True):
        r={'fm_Hz':int(fm),'n_sweeps':len(g)}
        for key in SIDEBANDS:
            mu=float(g[key].mean());sd=float(g[key].std(ddof=1))
            r.update({f'{key}_mean':mu,f'{key}_sample_SD':sd,
                      f'{key}_relative_SD_percent':100*sd/abs(mu) if mu else np.nan})
        for key in ('H_L','H_C','H_U'):
            r.update({f'{key}_mean':float(g[key].mean()),f'{key}_sample_SD':float(g[key].std(ddof=1))})
        for key in ('CH1_LSB_snr_dB','CH1_USB_snr_dB','CH2_LSB_snr_dB','CH2_USB_snr_dB'):
            r.update({f'{key}_min':float(g[key].min()),f'{key}_mean':float(g[key].mean())})
        rows.append(r)
    return pd.DataFrame(rows)


def process_formal(name, out):
    data=ROOT/'data'/name; output=out/name;output.mkdir(parents=True,exist_ok=True)
    raw=pd.read_csv(data/'formal_am_raw.csv')
    reference=pd.read_csv(data/'formal_am_summary.csv').sort_values('fm_Hz')
    ends=pd.read_csv(data/'formal_am_end_checks.csv')
    count,frequencies=FORMAL[name]
    if len(raw)!=count or raw.fm_Hz.nunique()!=frequencies or len(ends)!=5:
        raise ValueError(f'Unexpected row counts in {name}')
    if raw.duplicated(['sweep_number','fm_Hz']).any() or not (raw.groupby('fm_Hz').size()==5).all():
        raise ValueError(f'Missing or duplicate sweep records in {name}')
    for col in ('CH1_clipped','CH2_clipped'):
        if raw[col].any():raise ValueError(f'Clipped record in {name}')
    for line in ('LSB','carrier','USB'):
        if not raw[f'{line}_same_bin_CH1_CH2'].all():raise ValueError('Channel FFT bins disagree')
    assert_close(raw.T_lower,raw.H_L/raw.H_C,'lower normalization')
    assert_close(raw.T_upper,raw.H_U/raw.H_C,'upper normalization')
    assert_close(raw.T_mean,(raw.T_lower+raw.T_upper)/2,'mean sideband transfer')
    summary=formal_summary(raw)
    for key in reference.columns:
        assert_close(summary[key],reference[key],f'{name}/{key}')
    for key in SIDEBANDS:
        initial=ends[f'initial_{key}'];end=ends[f'endcheck_{key}']
        assert_close(ends[f'{key}_drift'],end-initial,f'{key} drift')
        assert_close(ends[f'{key}_drift_percent'],100*(end-initial)/initial,f'{key} percent drift')
        starts=raw[raw.fm_Hz==1000].set_index('sweep_number')[key]
        assert_close(initial,starts.loc[ends.sweep_number].to_numpy(),f'{key} end-check baseline')
    summary.to_csv(output/'formal_am_summary.csv',index=False)
    for filename,mode in [('T_lower_upper.png','sides'),('T_mean.png','mean'),('sweep_comparison.png','sweeps'),
                          ('repeatability.png','sd'),('endcheck_drift.png','end')]:
        fig,ax=plt.subplots(figsize=(9,5.4))
        if mode=='sides':
            for sw,g in raw.groupby('sweep_number'):
                line,=ax.plot(g.fm_Hz/1000,g.T_lower,'o-',ms=3,label=f'Sweep {sw} lower')
                ax.plot(g.fm_Hz/1000,g.T_upper,'s--',ms=3,color=line.get_color(),label=f'Sweep {sw} upper')
            ylabel='Normalized sideband transfer';title='Lower and upper sidebands'
        elif mode=='mean':
            ax.errorbar(summary.fm_Hz/1000,summary.T_mean_mean,yerr=summary.T_mean_sample_SD,fmt='o-',capsize=3)
            ylabel='T_mean';title='Across five sweeps: mean +/- 1 sample SD'
        elif mode=='sweeps':
            for sw,g in raw.groupby('sweep_number'):ax.plot(g.fm_Hz/1000,g.T_mean,'o-',ms=3,label=f'Sweep {sw}')
            ylabel='T_mean';title='Individual sweeps'
        elif mode=='sd':
            for k in SIDEBANDS:ax.plot(summary.fm_Hz/1000,summary[f'{k}_relative_SD_percent'],'o-',ms=3,label=k)
            ylabel='Relative sample SD (%)';title='Between-sweep repeatability'
        else:
            for k in SIDEBANDS:ax.plot(ends.sweep_number,ends[f'{k}_drift_percent'],'o-',label=k)
            ax.axhline(0,color='black',lw=.7);ylabel='Drift from initial 1 kHz (%)';title='End checks'
        ax.set(xlabel='Sweep' if mode=='end' else 'Modulation frequency (kHz)',ylabel=ylabel,title=f'{name}: {title}')
        ax.grid(alpha=.25)
        if mode!='mean':ax.legend(fontsize=7,ncol=2)
        save(fig,output/filename)
    return {'formal_measurements':len(raw),'frequency_points':frequencies,'end_checks':len(ends),
            'minimum_sideband_SNR_dB':float(raw[['CH1_LSB_snr_dB','CH1_USB_snr_dB','CH2_LSB_snr_dB','CH2_USB_snr_dB']].min().min()),
            'max_T_mean_sample_SD':float(summary.T_mean_sample_SD.max()),
            'max_absolute_end_drift_percent':float(ends[[k+'_drift_percent' for k in SIDEBANDS]].abs().to_numpy().max())}


def process_drive(name, out):
    filename='drive_dependence_support.csv' if name=='drive_support' else 'drive_dependence_results.csv'
    frame=pd.read_csv(ROOT/'data'/name/filename)
    if not frame.measurement_valid.all():raise ValueError('Invalid drive characterization row')
    expected=18 if name=='drive_support' else 21
    if len(frame)!=expected or frame.duplicated(['amplitude_programmed_Vpp','repeat']).any():raise ValueError('Invalid drive record count')
    output=out/name;output.mkdir(parents=True,exist_ok=True)
    keys=('f0_measured_Hz','Hmax','BW_measured_Hz','Q_measured','Vin_resonance_Vpp','Vout_resonance_Vpk',
          'I_resonance_pk_A_from_Vout_over_R')
    result=[]
    for amp,g in frame.groupby('amplitude_programmed_Vpp',sort=True):
        if len(g)!=3:raise ValueError('Expected three repeats per drive')
        assert_close(g.I_resonance_pk_A_from_Vout_over_R,g.Vout_resonance_Vpk/5.1,'current from nominal R')
        r={'amplitude_programmed_Vpp':amp,'n_repeats':len(g),'n_reliable_complex_fits':int(g.complex_fit_reliable.sum())}
        for k in keys:r.update({k+'_mean':g[k].mean(),k+'_sample_SD':g[k].std(ddof=1)})
        # Unreliable low-drive fits remain in source tables but are excluded from fit summaries.
        valid=g[g.complex_fit_reliable]
        for k in ('Rs_fit_ohm','L_eff_fit_uH','C_eff_fit_nF'):
            r[k+'_reliable_mean']=valid[k].mean() if len(valid) else np.nan
        result.append(r)
    summary=pd.DataFrame(result);summary.to_csv(output/'drive_summary.csv',index=False)
    for k,scale,label in [('f0_measured_Hz',1e-3,'Phase-zero resonance (kHz)'),('Hmax',1,'Peak transfer'),
                           ('BW_measured_Hz',1e-3,'Bandwidth (kHz)'),('Q_measured',1,'Q')]:
        fig,ax=plt.subplots(figsize=(8,5))
        regular=summary.amplitude_programmed_Vpp!=12
        ax.errorbar(summary.loc[regular,'amplitude_programmed_Vpp'],summary.loc[regular,k+'_mean']*scale,
                    yerr=summary.loc[regular,k+'_sample_SD']*scale,fmt='o-',capsize=3,label='Mean +/- 1 sample SD, n=3')
        extreme=summary[~regular]
        if len(extreme):
            ax.errorbar(extreme.amplitude_programmed_Vpp,extreme[k+'_mean']*scale,yerr=extreme[k+'_sample_SD']*scale,
                        fmt='o',color='red',capsize=3,label='12 Vpp extreme test')
        ax.set(xlabel='Generator setting (Vpp)',ylabel=label,title=name);ax.grid(alpha=.25);ax.legend()
        save(fig,output/(k+'.png'))
    return {'valid_repeats':len(frame),'drive_settings':len(summary),'unreliable_low_drive_fits':int((~frame.complex_fit_reliable).sum())}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset',choices=['all',*FORMAL,'drive_support','drive_extreme'],default='all')
    parser.add_argument('--output',type=Path,default=ROOT/'outputs')
    args=parser.parse_args()
    out=args.output.resolve()
    # Published measurements and preserved figures are immutable inputs.
    protected=[ROOT/'data',ROOT/'figures',ROOT/'docs',ROOT/'acquisition',ROOT/'analysis',ROOT/'tools',ROOT/'.git']
    if out==ROOT or any(out.is_relative_to(p) or p.is_relative_to(out) for p in protected):
        parser.error('Choose an output folder separate from published inputs and source code.')
    names=[*FORMAL,'drive_support','drive_extreme'] if args.dataset=='all' else [args.dataset]
    report={name:(process_formal(name,out) if name in FORMAL else process_drive(name,out)) for name in names}
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
