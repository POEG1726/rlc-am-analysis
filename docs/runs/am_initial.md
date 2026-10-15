# Formal AM data collection summary

- Collection completed: 2026-09-26 (UTC date; time withheld)
- Formal acquisitions: **75 / 75** (five sweeps × 15 frequencies).
- End checks: **5 / 5**; accepted acquisitions total: **80 / 80**.
- Rejected attempts preserved: **0**.
- Formal acquisition time span: 2026-09-26 (UTC date; time withheld) to 2026-09-26 (UTC date; time withheld).

## Acquisition validity

- Sideband SNR across CH1/CH2 and LSB/USB: **85.1 to 96.0 dB** (local per-bin FFT estimate).
- ADC clipping in accepted records: **none**.
- FFT-bin consistency: all accepted measurements passed the CH1/CH2 same-bin check; waveform data are native 5,000,000-point, 50 MSa/s captures.
- No interpolation or acquisition-time averaging was applied. The summary uses the five separate sweeps.

## Agreement and repeatability

- Largest sample SD of T_mean across modulation frequencies: **0.000658589** at **1 kHz**.
- Largest absolute initial-to-end 1 kHz drift among T_lower, T_upper, and T_mean: **0.1492%**. See `processed/formal_am_end_checks.csv`.

| fm (kHz) | mean T_lower ± SD | mean T_upper ± SD | mean T_mean ± SD | max relative SD T_mean (%) |
|---:|---:|---:|---:|---:|
| 1 | 1.000366 ± 0.00066 | 0.998879 ± 0.00066 | 0.999623 ± 0.00066 | 0.06588 |
| 2 | 0.998702 ± 0.00041 | 0.995935 ± 0.00039 | 0.997319 ± 0.0004 | 0.03989 |
| 3 | 0.995354 ± 0.00022 | 0.991405 ± 0.00022 | 0.993379 ± 0.00022 | 0.02215 |
| 5 | 0.984665 ± 0.00051 | 0.978755 ± 0.0005 | 0.981710 ± 0.0005 | 0.05115 |
| 8 | 0.958633 ± 0.00041 | 0.950980 ± 0.00038 | 0.954806 ± 0.0004 | 0.04149 |
| 10 | 0.935380 ± 0.00047 | 0.927442 ± 0.00049 | 0.931411 ± 0.00048 | 0.05148 |
| 12 | 0.909214 ± 0.00022 | 0.901764 ± 0.00022 | 0.905489 ± 0.00022 | 0.02442 |
| 15 | 0.865087 ± 0.00023 | 0.859545 ± 0.00021 | 0.862316 ± 0.00022 | 0.02533 |
| 18 | 0.818334 ± 0.00019 | 0.815732 ± 0.00021 | 0.817033 ± 0.0002 | 0.02488 |
| 20 | 0.786755 ± 0.00028 | 0.786503 ± 0.00028 | 0.786629 ± 0.00028 | 0.03581 |
| 22 | 0.755486 ± 0.00027 | 0.757757 ± 0.00027 | 0.756621 ± 0.00027 | 0.03537 |
| 25 | 0.709848 ± 0.00027 | 0.716036 ± 0.00028 | 0.712942 ± 0.00028 | 0.03877 |
| 30 | 0.639619 ± 0.00029 | 0.652117 ± 0.0003 | 0.645868 ± 0.0003 | 0.04627 |
| 35 | 0.577440 ± 0.00014 | 0.595599 ± 0.00014 | 0.586520 ± 0.00014 | 0.02363 |
| 40 | 0.523660 ± 0.0002 | 0.546661 ± 0.00021 | 0.535161 ± 0.0002 | 0.03817 |

## 1 kHz end checks

| Sweep | T_lower drift | T_upper drift | T_mean drift |
|---:|---:|---:|---:|
| 1 | -0.000855799 (-0.0855%) | -0.00090484 (-0.09053%) | -0.000880319 (-0.08801%) |
| 2 | +0.000462247 (+0.04622%) | +0.000492016 (+0.04927%) | +0.000477131 (+0.04775%) |
| 3 | -0.00149391 (-0.1492%) | -0.00144899 (-0.145%) | -0.00147145 (-0.1471%) |
| 4 | +0.000806669 (+0.0807%) | +0.000906002 (+0.09078%) | +0.000856335 (+0.08573%) |
| 5 | -0.000218645 (-0.02186%) | -0.000194722 (-0.0195%) | -0.000206683 (-0.02068%) |

## Anomalies and suitability

- No point exceeded three times its across-sweep sample SD from that frequency mean.
- Suitability: **suitable for IA analysis**: all five sweeps and end checks are complete, no clipping was detected, and the minimum measured sideband SNR was 85.1 dB. Interpret the reported SDs as between-sweep repeatability; sideband SNR is estimated locally from FFT bins.

## Outputs

- `processed/formal_am_raw.csv` — one row per formal acquisition.
- `processed/formal_am_summary.csv` — per-frequency means and sample SDs across five sweeps.
- `processed/formal_am_end_checks.csv` — within-sweep initial-to-end drift.
- Native waveform archives are retained privately and are not included in this repository.
- `figures/` — requested plots.
