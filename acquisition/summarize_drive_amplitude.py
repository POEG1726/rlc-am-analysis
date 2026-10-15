from __future__ import annotations

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import characterize_rlc as base
ROOT=_bench.workspace
R_EXT=5.1
cal=pd.read_csv(ROOT/'scope_vertical_range_calibration.csv')

def range_factors(ch,signal_col):
 d=cal[cal.tested_channel==ch]
 pairs=[]
 for amp,g in d.groupby('amplitude_programmed_Vpp'):
  gr=[]
  for scale,gg in g.groupby('tested_scale_V_div'):
   values=gg[signal_col].to_numpy(float);gr.append((float(scale),float(values.mean()),float(values.std(ddof=1)/math.sqrt(len(values)))) )
  gr.sort()
  if len(gr)!=2:continue
  (lo,ml,sel),(hi,mh,seh)=gr
  ratio=mh/ml
  relse=math.sqrt((seh/mh)**2+(sel/ml)**2)
  pairs.append({'amplitude_Vpp':float(amp),'low_scale_V_div':lo,'high_scale_V_div':hi,'low_mean_Vpk':ml,'high_mean_Vpk':mh,'low_over_high_reading':ml/mh,'correction_low_relative_to_high':ratio,'relative_standard_error':relse})
 factors={round(float(s),6):1.0 for s in d.tested_scale_V_div.unique()}
 reference=max(factors)
 # Build relative factors from the largest tested range downward through adjacent ranges.
 factors={reference:1.0}
 for pair in sorted(pairs,key=lambda x:x['high_scale_V_div'],reverse=True):
  lo=round(pair['low_scale_V_div'],6);hi=round(pair['high_scale_V_div'],6)
  if hi in factors:factors[lo]=pair['correction_low_relative_to_high']*factors[hi]
 return factors,pairs

cf1,pair1=range_factors(1,'vin_Vpk');cf2,pair2=range_factors(2,'vout_Vpk')
calibration={'reference_scale_V_div':{'CH1':max(cf1),'CH2':max(cf2)},'correction_factors_to_reference':{'CH1':cf1,'CH2':cf2},'pair_measurements':{'CH1':pair1,'CH2':pair2},'method':'For each adjacent scale pair, hold the same generator amplitude and frequency, change only the tested channel scale, discard one post-change capture, and use the ratio of mean measured amplitudes (three repeats each) to map readings to the 200 mV/div reference range. This corrects relative range gain, not absolute oscilloscope calibration.'}
(ROOT/'vertical_scale_corrections.json').write_text(json.dumps(calibration,indent=2),encoding='utf-8')
def factor(ch,scale):
 keys=cf1 if ch==1 else cf2
 k=round(float(scale),6)
 if k not in keys:raise KeyError(f'No vertical range calibration for CH{ch} at {scale:g} V/div')
 return keys[k]

raw=pd.read_csv(ROOT/'drive_amplitude_frequency_scan.csv')
for ch in (1,2):raw[f'CH{ch}_scale_correction_factor']=raw[f'CH{ch}_vdiv_V_div'].map(lambda x:factor(ch,x))
raw['Vin_corrected_Vpk']=raw.vin_Vpk*raw.CH1_scale_correction_factor
raw['Vout_corrected_Vpk']=raw.vout_Vpk*raw.CH2_scale_correction_factor
hcor=raw.CH2_scale_correction_factor/raw.CH1_scale_correction_factor
raw['H_real_corrected']=raw.H_real*hcor;raw['H_imag_corrected']=raw.H_imag*hcor;raw['H_mag_corrected']=raw.H_mag*hcor
raw.to_csv(ROOT/'drive_amplitude_frequency_scan.csv',index=False)
summary=pd.read_csv(ROOT/'drive_amplitude_sweep.csv')
# Retain raw values before replacing the primary voltage/gain/current fields with range-corrected values.
for name in ['Vin_resonance_Vpk','Vout_resonance_Vpk','Vin_resonance_Vpp','Hmax','I_resonance_pk_A','I_resonance_rms_A']:
 if f'{name}_uncorrected' not in summary:summary[f'{name}_uncorrected']=summary[name]
summary['CH1_scale_correction_factor']=summary.CH1_scale_V_div.map(lambda x:factor(1,x))
summary['CH2_scale_correction_factor']=summary.CH2_scale_V_div.map(lambda x:factor(2,x))
summary['Vin_resonance_Vpk']=summary.Vin_resonance_Vpk_uncorrected*summary.CH1_scale_correction_factor
summary['Vout_resonance_Vpk']=summary.Vout_resonance_Vpk_uncorrected*summary.CH2_scale_correction_factor
summary['Vin_resonance_Vpp']=2*summary.Vin_resonance_Vpk
ratio=summary.CH2_scale_correction_factor/summary.CH1_scale_correction_factor
summary['Hmax']=summary.Hmax_uncorrected*ratio
summary['I_resonance_pk_A']=summary.Vout_resonance_Vpk/R_EXT
summary['I_resonance_rms_A']=summary.I_resonance_pk_A/np.sqrt(2)
summary['Rs_from_Hmax_ohm']=R_EXT/summary.Hmax
# Refit complex transfer with measured vertical-range correction applied point-by-point.
fit_rows=[]
for (amp,direction),g in raw.groupby(['amplitude_programmed_Vpp','direction']):
 try:
  fit_data=g.copy();fit_data['H_real']=fit_data['H_real_corrected'];fit_data['H_imag']=fit_data['H_imag_corrected'];fit_data['H_mag']=fit_data['H_mag_corrected'];fit=base.fit_complex(fit_data)
  rs,L,C=fit['Rs_fit_ohm'],fit['L_fit_H'],fit['C_fit_F'];f0fit=1/(2*np.pi*np.sqrt(L*C));bwfit=rs/(2*np.pi*L)
  hmean=float(g.H_mag_corrected.mean());relrmse=fit['complex_H_RMSE']/hmean
  ok=bool(fit['fit_success'] and relrmse<=0.02 and abs(f0fit-float(summary[(summary.amplitude_programmed_Vpp==amp)&(summary.direction==direction)].f0_measured_Hz.iloc[0]))<0.01*float(summary[(summary.amplitude_programmed_Vpp==amp)&(summary.direction==direction)].f0_measured_Hz.iloc[0]))
  row={'amplitude_programmed_Vpp':amp,'direction':direction,'Rs_fit_ohm':rs,'L_eff_fit_H':L,'L_eff_fit_uH':L*1e6,'C_eff_fit_F':C,'C_eff_fit_nF':C*1e9,'f0_complex_fit_Hz':f0fit,'BW_complex_fit_Hz':bwfit,'Q_complex_fit':f0fit/bwfit,'complex_H_RMSE':fit['complex_H_RMSE'],'complex_H_relative_RMSE':relrmse,'channel_delay_fit_s':fit.get('channel_delay_fit_s',math.nan),'complex_fit_success':fit['fit_success'],'complex_fit_reliable':ok,'complex_fit_n':fit['fit_n']}
 except Exception as e:
  row={'amplitude_programmed_Vpp':amp,'direction':direction,'Rs_fit_ohm':math.nan,'L_eff_fit_H':math.nan,'L_eff_fit_uH':math.nan,'C_eff_fit_F':math.nan,'C_eff_fit_nF':math.nan,'f0_complex_fit_Hz':math.nan,'BW_complex_fit_Hz':math.nan,'Q_complex_fit':math.nan,'complex_H_RMSE':math.nan,'complex_H_relative_RMSE':math.nan,'channel_delay_fit_s':math.nan,'complex_fit_success':False,'complex_fit_reliable':False,'complex_fit_n':0,'fit_error':repr(e)}
 fit_rows.append(row)
fits=pd.DataFrame(fit_rows)
for col in ['Rs_fit_ohm','L_eff_fit_H','L_eff_fit_uH','C_eff_fit_F','C_eff_fit_nF','f0_complex_fit_Hz','BW_complex_fit_Hz','Q_complex_fit','complex_H_RMSE','complex_H_relative_RMSE','channel_delay_fit_s','complex_fit_success','complex_fit_reliable','complex_fit_n']:
 old=f'{col}_uncorrected'
 if old not in summary and col in summary:summary[old]=summary[col]
 summary=summary.drop(columns=[col],errors='ignore').merge(fits[['amplitude_programmed_Vpp','direction',col]],on=['amplitude_programmed_Vpp','direction'],how='left')
summary.to_csv(ROOT/'drive_amplitude_sweep.csv',index=False)
# Harmonic levels in dBc are scale invariant; also provide corrected absolute line amplitudes.
harm=pd.read_csv(ROOT/'drive_amplitude_harmonics.csv')
harm['CH1_scale_correction_factor']=harm.CH1_vdiv_V_div.map(lambda x:factor(1,x));harm['CH2_scale_correction_factor']=harm.CH2_vdiv_V_div.map(lambda x:factor(2,x))
for ch in (1,2):
 for nm in ['fundamental','H2','H3']:
  col=f'CH{ch}_{nm}_Vpk'
  harm[col+'_corrected']=harm[col]*harm[f'CH{ch}_scale_correction_factor']
harm.to_csv(ROOT/'drive_amplitude_harmonics.csv',index=False)
# Plots, ascending and descending kept as separate point series.
marker={'ascending':'o','descending':'s'};color={'ascending':'C0','descending':'C1'}
def plot_xy(xcol,ycol,xlabel,ylabel,title,path,y_scale=1):
 fig,ax=plt.subplots(figsize=(8.5,5.4))
 for dr,g in summary.groupby('direction'):
  g=g.sort_values('amplitude_programmed_Vpp');ax.plot(g[xcol],g[ycol]/y_scale,marker=marker[dr],ms=4,lw=1,label=dr,color=color[dr])
 ax.set(xlabel=xlabel,ylabel=ylabel,title=title);ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(ROOT/path,dpi=170);plt.close(fig)
plot_xy('amplitude_programmed_Vpp','f0_measured_Hz','Generator setting (Vpp)','Phase-zero resonance (kHz)','Resonance vs programmed drive','f0_vs_programmed_vpp.png',1000)
plot_xy('Vin_resonance_Vpp','f0_measured_Hz','Actual CH1 input at resonance (Vpp)','Phase-zero resonance (kHz)','Resonance vs measured input','f0_vs_actual_vin.png',1000)
plot_xy('I_resonance_pk_A','f0_measured_Hz','Resonance current (A peak)','Phase-zero resonance (kHz)','Resonance vs measured current','f0_vs_resonance_current.png',1000)
plot_xy('amplitude_programmed_Vpp','Hmax','Generator setting (Vpp)','Corrected Hmax','Peak transfer vs programmed drive','hmax_vs_programmed_vpp.png')
plot_xy('I_resonance_pk_A','Hmax','Resonance current (A peak)','Corrected Hmax','Peak transfer vs measured current','hmax_vs_current.png')
# Corrected fitted parameters with low-fit-quality points still visible but documented.
fig,axs=plt.subplots(3,1,figsize=(8.7,10),sharex=True)
for dr,g in summary.groupby('direction'):
 g=g.sort_values('amplitude_programmed_Vpp')
 axs[0].plot(g.amplitude_programmed_Vpp,g.Rs_fit_ohm,marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr])
 axs[1].plot(g.amplitude_programmed_Vpp,g.L_eff_fit_uH,marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr])
 axs[2].plot(g.amplitude_programmed_Vpp,g.C_eff_fit_nF,marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr])
axs[0].set_ylabel('Rs fit (Ω)');axs[1].set_ylabel('L_eff fit (µH)');axs[2].set_ylabel('C_eff fit (nF)');axs[2].set_xlabel('Generator setting (Vpp)')
for ax in axs:ax.grid(alpha=.25);ax.legend()
fig.suptitle('Range-corrected complex-fit parameters (hollow confidence not implied)');fig.tight_layout();fig.savefig(ROOT/'fitted_parameters_vs_amplitude.png',dpi=170);plt.close(fig)
fig,axs=plt.subplots(2,1,figsize=(8.5,7),sharex=True)
for dr,g in summary.groupby('direction'):
 g=g.sort_values('amplitude_programmed_Vpp');axs[0].plot(g.amplitude_programmed_Vpp,g.BW_measured_Hz/1000,marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr]);axs[1].plot(g.amplitude_programmed_Vpp,g.Q_measured,marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr])
axs[0].set_ylabel('BW (kHz)');axs[1].set_ylabel('Q');axs[1].set_xlabel('Generator setting (Vpp)')
for ax in axs:ax.grid(alpha=.25);ax.legend()
fig.suptitle('Bandwidth and Q vs drive amplitude');fig.tight_layout();fig.savefig(ROOT/'bandwidth_q_vs_amplitude.png',dpi=170);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(10,8),sharex=True)
for ax,ch,harmonic in [(axs[0,0],1,'H2'),(axs[0,1],1,'H3'),(axs[1,0],2,'H2'),(axs[1,1],2,'H3')]:
 col=f'CH{ch}_{harmonic}_dBc'
 for dr,g in summary.groupby('direction'):
  g=g.sort_values('amplitude_programmed_Vpp');ax.plot(g.amplitude_programmed_Vpp,g[col],marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr])
 ax.set_title(f'CH{ch} {harmonic}');ax.set_ylabel('Level (dBc)');ax.grid(alpha=.25);ax.legend()
for ax in axs[1,:]:ax.set_xlabel('Generator setting (Vpp)')
fig.suptitle('Harmonics at the phase-zero frequency');fig.tight_layout();fig.savefig(ROOT/'harmonics_vs_amplitude.png',dpi=170);plt.close(fig)
# Ascending minus descending differences at each matched amplitude.
a=summary[summary.direction=='ascending'].set_index('amplitude_programmed_Vpp').sort_index();d=summary[summary.direction=='descending'].set_index('amplitude_programmed_Vpp').sort_index()
comp=[('f0_measured_Hz','Δf0 (kHz)',1000),('Hmax','ΔHmax',1),('BW_measured_Hz','ΔBW (kHz)',1000),('Q_measured','ΔQ',1),('Rs_fit_ohm','ΔRs (Ω)',1),('L_eff_fit_uH','ΔL_eff (µH)',1),('C_eff_fit_nF','ΔC_eff (nF)',1)]
fig,axs=plt.subplots(len(comp),1,figsize=(8.5,14),sharex=True)
for ax,(k,label,scale) in zip(axs,comp):ax.plot(a.index,(a[k]-d[k])/scale,'o-',ms=3);ax.axhline(0,color='k',lw=.6);ax.set_ylabel(label);ax.grid(alpha=.25)
axs[-1].set_xlabel('Generator setting (Vpp)');fig.suptitle('Ascending minus descending');fig.tight_layout();fig.savefig(ROOT/'ascending_descending_comparison.png',dpi=170);plt.close(fig)
# Key statistics for the report.
mean=summary.groupby('amplitude_programmed_Vpp').mean(numeric_only=True)
low=mean.loc[[.1,.15,.2]];high=mean.loc[[.1,.2,.3,.4,.5,.6]];reliable=summary[summary.complex_fit_reliable]
maxdiff={k:float((a[k]-d[k]).abs().max()) for k,_,_ in comp}
print('FACTORS CH1',cf1);print('FACTORS CH2',cf2)
print('LOW 0.1-0.2 range:',{k:(float(low[k].min()),float(low[k].max())) for k in ['f0_measured_Hz','Hmax','BW_measured_Hz','Q_measured','Rs_fit_ohm','L_eff_fit_uH','C_eff_fit_nF']})
print('ALL low/high rows:');print(mean[['f0_measured_Hz','Hmax','BW_measured_Hz','Q_measured','Vin_resonance_Vpp','I_resonance_pk_A','Rs_fit_ohm','L_eff_fit_uH','C_eff_fit_nF','complex_H_relative_RMSE','complex_fit_reliable','CH1_H2_dBc','CH1_H3_dBc','CH2_H2_dBc','CH2_H3_dBc']].to_string())
print('RELIABLE amps', sorted(reliable.amplitude_programmed_Vpp.unique()));print('MAX up-down',maxdiff)



# Fitted parameter current plot.
fig,axs=plt.subplots(3,1,figsize=(8.7,10),sharex=True)
for dr,g in summary.groupby('direction'):
 g=g.sort_values('I_resonance_pk_A')
 axs[0].plot(g.I_resonance_pk_A*1000,g.Rs_fit_ohm,marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr])
 axs[1].plot(g.I_resonance_pk_A*1000,g.L_eff_fit_uH,marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr])
 axs[2].plot(g.I_resonance_pk_A*1000,g.C_eff_fit_nF,marker=marker[dr],ms=3,lw=1,label=dr,color=color[dr])
axs[0].set_ylabel('Rs fit (Ω)');axs[1].set_ylabel('L_eff fit (µH)');axs[2].set_ylabel('C_eff fit (nF)');axs[2].set_xlabel('Resonance current (mA peak)')
for ax in axs:ax.grid(alpha=.25);ax.legend()
fig.suptitle('Range-corrected complex-fit parameters vs resonance current');fig.tight_layout();fig.savefig(ROOT/'fitted_parameters_vs_current.png',dpi=170);plt.close(fig)
