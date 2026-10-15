from __future__ import annotations

# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import json, math, re, time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyvisa
import matplotlib.pyplot as plt

import characterize_rlc as base

ROOT=_bench.workspace
R_EXT=5.1
AMPLITUDES=np.array([0.10,0.15,0.20,0.30,0.40,0.50,0.60,0.75,1.00,1.25,1.50,1.75,2.00,2.50,3.00,3.50,4.00,4.50,5.00])
COARSE=np.arange(450_000.0,550_000.1,5_000.0)
NAMES={
 'summary':'drive_amplitude_sweep.csv','raw':'drive_amplitude_frequency_scan.csv','harmonics':'drive_amplitude_harmonics.csv','log':'drive_amplitude_log.txt'
}
SUPPORTED_SCALES=np.array([m*10.0**e for e in range(-3,1) for m in (1,2,5)],dtype=float)
RUN_LOG=[]

def log(s):
 print(s,flush=True);RUN_LOG.append(str(s))

def nearest_safe_scale(vpk,target_div=2.5):
 target=max(vpk/target_div,1e-3)
 choices=SUPPORTED_SCALES[SUPPORTED_SCALES>=target]
 return float(choices[0] if len(choices) else SUPPORTED_SCALES[-1])

def set_scales(scope,amp_vpp):
 # Pilot calibration: CH1 is about 0.063*Vpp peak and CH2 about 0.046*Vpp peak.
 s1=nearest_safe_scale(0.063*amp_vpp)
 s2=nearest_safe_scale(0.046*amp_vpp)
 scope.write(f':CHAN1:SCAL {s1:.9g}');scope.write(f':CHAN2:SCAL {s2:.9g}')
 actual1=float(scope.query(':CHAN1:SCAL?'));actual2=float(scope.query(':CHAN2:SCAL?'))
 return actual1,actual2

def auto_range_at(scope,gen,amp_vpp,fcenter):
 # Verify the estimate with actual waveform codes and alter one scale step if unsafe.
 s1,s2=set_scales(scope,amp_vpp); row=None
 for _ in range(3):
  row=base.measure(scope,gen,float(fcenter),'range_check')
  scales=[float(scope.query(':CHAN1:SCAL?')),float(scope.query(':CHAN2:SCAL?'))]
  occ=[max(abs(row['CH1_code_min']),abs(row['CH1_code_max']))/7680.0,
       max(abs(row['CH2_code_min']),abs(row['CH2_code_max']))/7680.0]
  changed=False
  for ch,(scale,div) in enumerate(zip(scales,occ),start=1):
   if div>3.55:
    candidates=SUPPORTED_SCALES[SUPPORTED_SCALES>scale]
    if len(candidates): scope.write(f':CHAN{ch}:SCAL {candidates[0]:.9g}');changed=True
   elif div<0.75:
    candidates=SUPPORTED_SCALES[SUPPORTED_SCALES<scale]
    if len(candidates): scope.write(f':CHAN{ch}:SCAL {candidates[-1]:.9g}');changed=True
  if not changed: break
 s1=float(scope.query(':CHAN1:SCAL?'));s2=float(scope.query(':CHAN2:SCAL?'))
 return row,s1,s2

def phase_zero(frame):
 d=frame.sort_values('frequency_measured_Hz').copy()
 f=d.frequency_measured_Hz.to_numpy(); ph=np.unwrap(np.deg2rad(d.H_phase_deg.to_numpy()))
 changes=np.flatnonzero(ph[:-1]*ph[1:]<=0)
 if not len(changes): raise RuntimeError('No phase-zero crossing in coarse sweep')
 i=int(changes[np.argmin(abs(f[changes]-495000))])
 lo=max(0,i-2);hi=min(len(d),i+4)
 ff=f[lo:hi];pp=np.rad2deg(ph[lo:hi]);center=float(np.mean(ff))
 slope,intercept=np.polyfit(ff-center,pp,1)
 fm=center-intercept/slope
 # Also translate local zero to programmed setting, minimizing generator clock-axis dependence.
 fs=d.frequency_set_Hz.to_numpy()[lo:hi];cs=float(np.mean(fs));slope_s,intercept_s=np.polyfit(fs-cs,pp,1);fc=cs-intercept_s/slope_s
 rms=float(np.sqrt(np.mean((pp-(slope*(ff-center)+intercept))**2)))
 return float(fm),float(fc),rms

def crossing_guess(frame,target,left=True):
 d=frame.sort_values('frequency_measured_Hz');f=d.frequency_measured_Hz.to_numpy();h=d.H_mag.to_numpy()
 ip=int(np.argmax(h)); choices=[]
 if left:
  for i in range(ip):
   if h[i]<target<=h[i+1]:choices.append(i)
 else:
  for i in range(ip,len(h)-1):
   if h[i]>=target>h[i+1]:choices.append(i)
 if not choices:return math.nan
 i=choices[-1] if left else choices[0]
 return float(np.interp(target,[h[i],h[i+1]],[f[i],f[i+1]]))

def direct_crossings(frame,hmax):
 d=frame.sort_values('frequency_measured_Hz'); f=d.frequency_measured_Hz.to_numpy();h=d.H_mag.to_numpy(); target=hmax/math.sqrt(2); ip=int(np.argmin(abs(h-hmax)))
 left=[i for i in range(ip) if h[i]<target<=h[i+1]]
 right=[i for i in range(ip,len(h)-1) if h[i]>=target>h[i+1]]
 if not left or not right:return math.nan,math.nan
 i=left[-1];fl=np.interp(target,[h[i],h[i+1]],[f[i],f[i+1]])
 i=right[0];fu=np.interp(target,[h[i+1],h[i]],[f[i+1],f[i]])
 return float(fl),float(fu)

def spectrum_peak(freq,amp,target,width):
 ids=np.flatnonzero((freq>=target-width)&(freq<=target+width))
 if not len(ids):return math.nan
 k=int(ids[np.argmax(amp[ids])]);return float(amp[k])

def measure_harmonics(scope,f0,on_capture=None):
 scope.write(':TRIG:RUN');time.sleep(.50);scope.write(':TRIG:STOP');time.sleep(.02)
 x1,dt1,m1=base.acquire_channel(scope,1);x2,dt2,m2=base.acquire_channel(scope,2)
 fs=1/dt1;n=min(len(x1),len(x2));freq=np.fft.rfftfreq(n,dt1);window=np.hanning(n);cg=window.mean()
 def amps(x):
  x=x[:n]-np.mean(x[:n]);a=2*np.abs(np.fft.rfft(x*window))/(n*cg);a[0]/=2;return a
 a1=amps(x1);a2=amps(x2);df=fs/n
 targets=[f0,2*f0,3*f0];width=max(1500.0,3*df)
 vals=[]
 for a in (a1,a2):vals.append([spectrum_peak(freq,a,t,width) for t in targets])
 def dbc(harm,fund):return 20*np.log10(max(harm,1e-15)/max(fund,1e-15))
 row={'harmonic_fundamental_Hz':f0,'frequency_resolution_Hz':df,
      'CH1_fundamental_Vpk':vals[0][0],'CH1_H2_Vpk':vals[0][1],'CH1_H3_Vpk':vals[0][2],
      'CH1_H2_dBc':dbc(vals[0][1],vals[0][0]),'CH1_H3_dBc':dbc(vals[0][2],vals[0][0]),
      'CH2_fundamental_Vpk':vals[1][0],'CH2_H2_Vpk':vals[1][1],'CH2_H3_Vpk':vals[1][2],
      'CH2_H2_dBc':dbc(vals[1][1],vals[1][0]),'CH2_H3_dBc':dbc(vals[1][2],vals[1][0]),
      'CH1_vdiv_V_div':m1['vdiv_V_div'],'CH2_vdiv_V_div':m2['vdiv_V_div'],
      'CH1_code_min':m1['code_min'],'CH1_code_max':m1['code_max'],'CH2_code_min':m2['code_min'],'CH2_code_max':m2['code_max']}
 return row

def one_sweep(scope,gen,amp,direction,order,center_guess,raw_path,on_capture=None,on_programmed=None,on_harmonic=None):
 gen.write('C2:MDWV STATE,OFF')
 gen.write(f'C2:BSWV WVTP,SINE,FRQ,{center_guess:.3f}HZ,AMP,{amp:.6g}V,OFST,0V')
 gen.write('C2:OUTP ON,LOAD,HZ');time.sleep(.30)
 readback=gen.query('C2:BSWV?').strip()
 m=re.search(r'AMP,([0-9.+\-Ee]+)V',readback);amp_set=float(m.group(1)) if m else amp
 if on_programmed is not None:on_programmed(amp_set,readback)
 pilot,s1,s2=auto_range_at(scope,gen,amp_set,center_guess)
 log(f'[{order:02d}/39] {direction} {amp_set:.2f} Vpp scales CH1={s1:g} V/div CH2={s2:g} V/div; range check Vin={pilot["vin_Vpk"]:.6g}, Vout={pilot["vout_Vpk"]:.6g}')
 rows=[]
 def add(freq,label):
  r=base.measure(scope,gen,float(freq),label);r.update({'amplitude_programmed_Vpp':amp_set,'direction':direction,'sweep_order':order});rows.append(r)
  if on_capture is not None:on_capture(r)
  return r
 for f in COARSE:add(f,'coarse')
 coarse=pd.DataFrame(rows)
 f0act,f0set,phase_rms=phase_zero(coarse)
 fine_set=np.arange(max(450000,round((f0set-2000)/250)*250),min(550000,f0set+2000)+125,250,dtype=float)
 for f in fine_set:add(f,'phase_refinement')
 data=pd.DataFrame(rows);fine=data[data.label=='phase_refinement']
 f0act,f0set,phase_rms=phase_zero(fine)
 # Refine both half-power crossings using the coarse response as a bracket.
 hmax=float(fine.H_mag.max());target=hmax/math.sqrt(2)
 gl=crossing_guess(coarse,target,True);gu=crossing_guess(coarse,target,False)
 if not np.isfinite(gl+gu):raise RuntimeError(f'Could not bracket -3dB points at {amp} Vpp')
 bw_freqs=[]
 for est in (gl,gu):
  bw_freqs.extend(np.arange(max(450000,round((est-2500)/1000)*1000),min(550000,est+2500)+500,1000,dtype=float))
 for f in sorted(set(bw_freqs)):add(f,'bandwidth_refinement')
 data=pd.DataFrame(rows);hmax=float(data[data.label=='phase_refinement'].H_mag.max())
 fl,fu=direct_crossings(data,hmax);bw=fu-fl;q=f0act/bw
 # Direct near-resonance channels and current at phase-zero programmed point.
 resonance=add(f0set,'resonance_measurement')
 current_pk=resonance['vout_Vpk']/R_EXT
 # Harmonics at resonance, with current auto-selected vertical scales retained.
 harm=measure_harmonics(scope,f0act);harm.update({'amplitude_programmed_Vpp':amp_set,'direction':direction,'sweep_order':order,'Vin_resonance_Vpk':resonance['vin_Vpk'],'Vout_resonance_Vpk':resonance['vout_Vpk'],'I_resonance_pk_A':current_pk})
 if on_harmonic is not None:on_harmonic(harm)
 broad=data[data.label.isin(['coarse','phase_refinement','bandwidth_refinement'])]
 try:
  fit=base.fit_complex(broad);rs=fit['Rs_fit_ohm'];L=fit['L_fit_H'];C=fit['C_fit_F'];f0fit=1/(2*np.pi*np.sqrt(L*C));bwfit=rs/(2*np.pi*L);qfit=f0fit/bwfit
  fit_ok=bool(fit['fit_success'] and fit['complex_H_RMSE']<0.02 and abs(f0fit-f0act)<0.02*f0act)
 except Exception as e:
  fit={'Rs_fit_ohm':math.nan,'L_fit_H':math.nan,'C_fit_F':math.nan,'complex_H_RMSE':math.nan,'fit_success':False,'fit_message':repr(e)};rs=L=C=f0fit=bwfit=qfit=math.nan;fit_ok=False
 summary={'amplitude_programmed_Vpp':amp_set,'direction':direction,'sweep_order':order,'generator_readback':readback,
  'CH1_scale_V_div':s1,'CH2_scale_V_div':s2,'f0_measured_Hz':f0act,'f0_programmed_Hz':f0set,'phase_fit_RMS_deg':phase_rms,
  'Hmax':hmax,'f_lower_3dB_Hz':fl,'f_upper_3dB_Hz':fu,'BW_measured_Hz':bw,'Q_measured':q,
  'Vin_resonance_Vpk':resonance['vin_Vpk'],'Vin_resonance_Vpp':2*resonance['vin_Vpk'],
  'Vout_resonance_Vpk':resonance['vout_Vpk'],'I_resonance_pk_A':current_pk,'I_resonance_rms_A':current_pk/math.sqrt(2),
  'Rs_fit_ohm':rs,'L_eff_fit_H':L,'L_eff_fit_uH':L*1e6,'C_eff_fit_F':C,'C_eff_fit_nF':C*1e9,
  'f0_complex_fit_Hz':f0fit,'BW_complex_fit_Hz':bwfit,'Q_complex_fit':qfit,'complex_H_RMSE':fit.get('complex_H_RMSE',math.nan),
  'complex_fit_success':fit.get('fit_success',False),'complex_fit_reliable':fit_ok,'complex_fit_n':fit.get('fit_n',0),
  'CH1_H2_dBc':harm['CH1_H2_dBc'],'CH1_H3_dBc':harm['CH1_H3_dBc'],'CH2_H2_dBc':harm['CH2_H2_dBc'],'CH2_H3_dBc':harm['CH2_H3_dBc'],
  'CH1_peak_divisions':max(abs(resonance['CH1_code_min']),abs(resonance['CH1_code_max']))/7680,
  'CH2_peak_divisions':max(abs(resonance['CH2_code_min']),abs(resonance['CH2_code_max']))/7680,
  'start_utc':datetime.now(timezone.utc).isoformat()}
 for r in rows:r['amplitude_programmed_Vpp']=amp_set;r['direction']=direction;r['sweep_order']=order
 if on_capture is None:
  append=Path(raw_path).exists();pd.DataFrame(rows).to_csv(raw_path,mode='a' if append else 'w',header=not append,index=False)
 return summary,harm,rows

def plot_summary(df,downs):
 # Collapse ascending/descending series into repeat-style plot lines; retain points separately.
 order={'ascending':'o','descending':'s'}; colors={'ascending':'C0','descending':'C1'}
 def plot_xy(xkey,ykey,xlab,ylab,title,file,logx=False):
  fig,ax=plt.subplots(figsize=(8.5,5.4))
  for direction,g in df.groupby('direction'):
   g=g.sort_values('amplitude_programmed_Vpp');ax.plot(g[xkey],g[ykey],marker=order[direction],ms=4,lw=1.1,color=colors[direction],label=direction)
  ax.set(xlabel=xlab,ylabel=ylab,title=title);ax.grid(alpha=.25);ax.legend();
  if logx:ax.set_xscale('log')
  fig.tight_layout();fig.savefig(ROOT/file,dpi=170);plt.close(fig)
 # Resonance against programmed amplitude.
 fig,ax=plt.subplots(figsize=(8.5,5.4))
 for direction,g in df.groupby('direction'):
  g=g.sort_values('amplitude_programmed_Vpp');ax.plot(g.amplitude_programmed_Vpp,g.f0_measured_Hz/1000,marker=order[direction],ms=4,lw=1.1,color=colors[direction],label=direction)
 ax.set(xlabel='Generator setting (Vpp)',ylabel='Phase-zero resonance (kHz)',title='Resonance vs programmed drive');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(ROOT/'f0_vs_programmed_vpp.png',dpi=170);plt.close(fig)
 fig,ax=plt.subplots(figsize=(8.5,5.4))
 for direction,g in df.groupby('direction'):
  ax.plot(g.Vin_resonance_Vpp,g.f0_measured_Hz/1000,marker=order[direction],ms=4,lw=1.1,color=colors[direction],label=direction)
 ax.set(xlabel='Actual CH1 at resonance (Vpp)',ylabel='Phase-zero resonance (kHz)',title='Resonance vs measured input voltage');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(ROOT/'f0_vs_actual_vin.png',dpi=170);plt.close(fig)
 fig,ax=plt.subplots(figsize=(8.5,5.4))
 for direction,g in df.groupby('direction'):
  ax.plot(g.I_resonance_pk_A*1000,g.f0_measured_Hz/1000,marker=order[direction],ms=4,lw=1.1,color=colors[direction],label=direction)
 ax.set(xlabel='Resonance current (mA peak)',ylabel='Phase-zero resonance (kHz)',title='Resonance vs measured current');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(ROOT/'f0_vs_resonance_current.png',dpi=170);plt.close(fig)
 fig,ax=plt.subplots(figsize=(8.5,5.4))
 for direction,g in df.groupby('direction'):
  ax.plot(g.I_resonance_pk_A*1000,g.Hmax,marker=order[direction],ms=4,lw=1.1,color=colors[direction],label=direction)
 ax.set(xlabel='Resonance current (mA peak)',ylabel='Hmax=|Vout/Vin|',title='Peak transfer vs resonance current');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(ROOT/'hmax_vs_current.png',dpi=170);plt.close(fig)
 fig,ax=plt.subplots(figsize=(8.5,5.4))
 for direction,g in df.groupby('direction'):
  g=g.sort_values('amplitude_programmed_Vpp');ax.plot(g.amplitude_programmed_Vpp,g.Hmax,marker=order[direction],ms=4,lw=1.1,color=colors[direction],label=direction)
 ax.set(xlabel='Generator setting (Vpp)',ylabel='Hmax=|Vout/Vin|',title='Peak transfer vs programmed drive');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(ROOT/'hmax_vs_programmed_vpp.png',dpi=170);plt.close(fig)
 fig,axs=plt.subplots(3,1,figsize=(8.5,10),sharex=True)
 for direction,g in df.groupby('direction'):
  g=g.sort_values('amplitude_programmed_Vpp')
  axs[0].plot(g.amplitude_programmed_Vpp,g.Rs_fit_ohm,marker=order[direction],ms=3,lw=1,label=direction,color=colors[direction])
  axs[1].plot(g.amplitude_programmed_Vpp,g.L_eff_fit_uH,marker=order[direction],ms=3,lw=1,color=colors[direction])
  axs[2].plot(g.amplitude_programmed_Vpp,g.C_eff_fit_nF,marker=order[direction],ms=3,lw=1,color=colors[direction])
 axs[0].set_ylabel('Rs fit (Ω)');axs[1].set_ylabel('L_eff fit (µH)');axs[2].set_ylabel('C_eff fit (nF)');axs[2].set_xlabel('Generator setting (Vpp)')
 for ax in axs:ax.grid(alpha=.25);ax.legend()
 fig.suptitle('Complex-fit parameters vs drive amplitude');fig.tight_layout();fig.savefig(ROOT/'fitted_parameters_vs_amplitude.png',dpi=170);plt.close(fig)
 fig,axs=plt.subplots(2,1,figsize=(8.5,7),sharex=True)
 for direction,g in df.groupby('direction'):
  g=g.sort_values('amplitude_programmed_Vpp');axs[0].plot(g.amplitude_programmed_Vpp,g.BW_measured_Hz/1000,marker=order[direction],ms=3,lw=1,label=direction,color=colors[direction]);axs[1].plot(g.amplitude_programmed_Vpp,g.Q_measured,marker=order[direction],ms=3,lw=1,label=direction,color=colors[direction])
 axs[0].set_ylabel('BW (kHz)');axs[1].set_ylabel('Q');axs[1].set_xlabel('Generator setting (Vpp)')
 for ax in axs:ax.grid(alpha=.25);ax.legend()
 fig.suptitle('Bandwidth and Q vs drive amplitude');fig.tight_layout();fig.savefig(ROOT/'bandwidth_q_vs_amplitude.png',dpi=170);plt.close(fig)
 fig,axs=plt.subplots(2,1,figsize=(8.5,7),sharex=True)
 for direction,g in df.groupby('direction'):
  g=g.sort_values('amplitude_programmed_Vpp')
  for ax,key,label in [(axs[0],'CH1_H2_dBc','CH1 2nd'),(axs[0],'CH1_H3_dBc','CH1 3rd'),(axs[1],'CH2_H2_dBc','CH2 2nd'),(axs[1],'CH2_H3_dBc','CH2 3rd')]:ax.plot(g.amplitude_programmed_Vpp,g[key],marker=order[direction],ms=3,lw=1,label=f'{direction} {label}',color=colors[direction] if '2nd' in label else ('C2' if direction=='ascending' else 'C3'))
 axs[0].set_ylabel('CH1 harmonic level (dBc)');axs[1].set_ylabel('CH2 harmonic level (dBc)');axs[1].set_xlabel('Generator setting (Vpp)')
 for ax in axs:ax.grid(alpha=.25);ax.legend(fontsize=8,ncol=2)
 fig.suptitle('2nd and 3rd harmonic levels at resonance');fig.tight_layout();fig.savefig(ROOT/'harmonics_vs_amplitude.png',dpi=170);plt.close(fig)
 # Ascending vs descending absolute difference across outputs
 asc=df[df.direction=='ascending'].set_index('amplitude_programmed_Vpp');desc=df[df.direction=='descending'].set_index('amplitude_programmed_Vpp')
 keys=[('f0_measured_Hz','Δf0 (kHz)',1000),('Hmax','ΔHmax',1),('BW_measured_Hz','ΔBW (kHz)',1000),('Q_measured','ΔQ',1),('Rs_fit_ohm','ΔRs (Ω)',1),('L_eff_fit_uH','ΔL_eff (µH)',1),('C_eff_fit_nF','ΔC_eff (nF)',1)]
 fig,axs=plt.subplots(len(keys),1,figsize=(8.5,14),sharex=True)
 for ax,(key,label,scale) in zip(axs,keys):ax.plot(asc.index,(asc[key]-desc[key])/scale,'o-',ms=3);ax.axhline(0,color='k',lw=.6);ax.set_ylabel(label);ax.grid(alpha=.25)
 axs[-1].set_xlabel('Generator setting (Vpp)');fig.suptitle('Ascending minus descending sweep');fig.tight_layout();fig.savefig(ROOT/'ascending_descending_comparison.png',dpi=170);plt.close(fig)

def main():
 rm=pyvisa.ResourceManager();scope=rm.open_resource(base.SCOPE_ADDR);gen=rm.open_resource(base.GEN_ADDR)
 scope.timeout=120000;gen.timeout=8000;scope.chunk_size=4*1024*1024
 summary_path=Path(NAMES['summary']);harm_path=Path(NAMES['harmonics']);raw_path=Path(NAMES['raw']);log_path=Path(NAMES['log'])
 summary_rows=pd.read_csv(summary_path).to_dict('records') if summary_path.exists() else []
 harm_rows=pd.read_csv(harm_path).to_dict('records') if harm_path.exists() else []
 completed={(str(r['direction']),round(float(r['amplitude_programmed_Vpp']),4)) for r in summary_rows}
 center=float(summary_rows[-1]['f0_programmed_Hz']) if summary_rows else 497150.0
 order=max([int(r['sweep_order']) for r in summary_rows],default=0)
 if log_path.exists():RUN_LOG.extend(log_path.read_text(encoding='utf-8').splitlines())
 try:
  scope.write(':TIM:SCAL 2E-4')
  log(f'Start UTC {datetime.now(timezone.utc).isoformat()}')
  log('Scope '+scope.query('*IDN?').strip());log('Generator '+gen.query('*IDN?').strip())
  log('Generator initial '+gen.query('C2:BSWV?').strip()+'; '+gen.query('C2:MDWV?').strip()+'; '+gen.query('C2:OUTP?').strip())
  log(f'Scope points={scope.query(":ACQ:POIN?").strip()} rate={scope.query(":ACQ:SRAT?").strip()}')
  if 'STATE,OFF' not in gen.query('C2:MDWV?').upper():gen.write('C2:MDWV STATE,OFF')
  sequence=[('ascending',float(a)) for a in AMPLITUDES]+[('descending',float(a)) for a in AMPLITUDES[::-1]]
  for direction,amp in sequence:
   if (direction,round(amp,4)) in completed:continue
   order+=1;time.sleep(.35)
   row,harm,raw=one_sweep(scope,gen,amp,direction,order,center,raw_path)
   summary_rows.append(row);harm_rows.append(harm);center=row['f0_programmed_Hz']
   pd.DataFrame(summary_rows).to_csv(summary_path,index=False);pd.DataFrame(harm_rows).to_csv(harm_path,index=False)
   log(f'complete {direction} {amp:.2f}Vpp: f0={row["f0_measured_Hz"]:.1f}Hz, Hmax={row["Hmax"]:.5f}, BW={row["BW_measured_Hz"]:.1f}Hz, Q={row["Q_measured"]:.3f}, Ipk={row["I_resonance_pk_A"]*1000:.2f}mA, Rs/L/C={row["Rs_fit_ohm"]:.4f}/{row["L_eff_fit_uH"]:.4f}/{row["C_eff_fit_nF"]:.4f}, harmonics CH1={row["CH1_H2_dBc"]:.1f}/{row["CH1_H3_dBc"]:.1f}dBc')
  summary=pd.DataFrame(summary_rows);summary.to_csv(summary_path,index=False);pd.DataFrame(harm_rows).to_csv(harm_path,index=False)
  plot_summary(summary,None)
  log('Sweep complete; wrote summary and figures.')
 finally:
  try:
   # Restore long record, original channel ranges, sine on CH2 at measured high-drive resonance, AM OFF.
   gen.write('C2:MDWV STATE,OFF');gen.write('C2:BSWV WVTP,SINE,FRQ,493900HZ,AMP,5V,OFST,0V');gen.write('C2:OUTP ON,LOAD,HZ')
   scope.write(':CHAN1:SCAL 2E-1');scope.write(':CHAN2:SCAL 1E-1');scope.write(':TIM:SCAL 5E-2');scope.write(':TRIG:RUN')
   log('Final generator '+gen.query('C2:BSWV?').strip()+'; '+gen.query('C2:MDWV?').strip()+'; '+gen.query('C2:OUTP?').strip())
   log(f'Final scope points={scope.query(":ACQ:POIN?").strip()} sample_rate={scope.query(":ACQ:SRAT?").strip()} timebase={scope.query(":TIM:SCAL?").strip()} CH1={scope.query(":CHAN1:IMP?").strip()} CH2={scope.query(":CHAN2:IMP?").strip()}')
   scope.close();gen.close();rm.close()
  except Exception as e:log('Restore error: '+repr(e))
  Path(NAMES['log']).write_text('\n'.join(RUN_LOG)+'\n',encoding='utf-8')

if __name__=='__main__':main()
