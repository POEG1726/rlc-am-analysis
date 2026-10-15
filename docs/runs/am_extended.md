# Formal AM data collection summary

- Collection completed: 2026-09-30 (UTC date; time withheld)
- Formal acquisitions: **155 / 155** (five sweeps × 31 frequencies).
- End checks: **5 / 5**; accepted acquisitions total: **160 / 160**.
- Rejected attempts preserved: **0**.
- Formal acquisition time span: 2026-09-30 (UTC date; time withheld) to 2026-09-30 (UTC date; time withheld).

## Acquisition validity

- Sideband SNR across CH1/CH2 and LSB/USB: **85.5 to 102.7 dB** (local per-bin FFT estimate).
- ADC clipping in accepted records: **none**.
- FFT-bin consistency: all accepted measurements passed the CH1/CH2 same-bin check; waveform data are native 5,000,000-point, 50 MSa/s captures.
- No interpolation or acquisition-time averaging was applied. The summary uses the five separate sweeps.

## Agreement and repeatability

- Largest sample SD of T_mean across modulation frequencies: **0.000510655** at **8 kHz**.
- Largest absolute initial-to-end 1 kHz drift among T_lower, T_upper, and T_mean: **0.07742%**. See `processed/formal_am_end_checks.csv`.

| fm (kHz) | mean T_lower ± SD | mean T_upper ± SD | mean T_mean ± SD | max relative SD T_mean (%) |
|---:|---:|---:|---:|---:|
| 1 | 0.996332 ± 0.00049 | 0.994198 ± 0.00046 | 0.995265 ± 0.00047 | 0.04746 |
| 2 | 0.994995 ± 0.00016 | 0.990980 ± 0.00017 | 0.992988 ± 0.00016 | 0.01645 |
| 3 | 0.992209 ± 0.00027 | 0.986404 ± 0.00028 | 0.989307 ± 0.00027 | 0.02769 |
| 5 | 0.981930 ± 0.00039 | 0.972995 ± 0.00042 | 0.977462 ± 0.0004 | 0.04123 |
| 8 | 0.956760 ± 0.00052 | 0.944625 ± 0.0005 | 0.950693 ± 0.00051 | 0.05371 |
| 10 | 0.933513 ± 0.00035 | 0.920424 ± 0.00032 | 0.926968 ± 0.00033 | 0.03572 |
| 12 | 0.907794 ± 0.00032 | 0.894646 ± 0.00034 | 0.901220 ± 0.00033 | 0.0369 |
| 15 | 0.863901 ± 0.00025 | 0.852215 ± 0.00023 | 0.858058 ± 0.00024 | 0.02809 |
| 18 | 0.817729 ± 0.00019 | 0.808853 ± 0.00024 | 0.813291 ± 0.00021 | 0.02618 |
| 20 | 0.786299 ± 0.00021 | 0.779830 ± 0.00026 | 0.783065 ± 0.00024 | 0.03017 |
| 22 | 0.755103 ± 0.00023 | 0.751280 ± 0.00025 | 0.753191 ± 0.00024 | 0.03195 |
| 25 | 0.709933 ± 0.00026 | 0.710288 ± 0.00025 | 0.710111 ± 0.00025 | 0.0359 |
| 30 | 0.640079 ± 0.00016 | 0.647364 ± 0.00014 | 0.643721 ± 0.00015 | 0.02304 |
| 35 | 0.578113 ± 0.00026 | 0.591690 ± 0.00031 | 0.584902 ± 0.00028 | 0.04856 |
| 40 | 0.523932 ± 0.00025 | 0.542943 ± 0.00028 | 0.533437 ± 0.00026 | 0.04908 |
| 45 | 0.476947 ± 0.00017 | 0.500508 ± 0.0002 | 0.488728 ± 0.00018 | 0.03692 |
| 50 | 0.436376 ± 9.7e-05 | 0.463759 ± 0.00011 | 0.450068 ± 0.0001 | 0.02308 |
| 60 | 0.370245 ± 0.00012 | 0.403393 ± 0.00016 | 0.386819 ± 0.00014 | 0.03601 |
| 70 | 0.319371 ± 8.2e-05 | 0.356583 ± 8.5e-05 | 0.337977 ± 8.4e-05 | 0.02472 |
| 80 | 0.279171 ± 9.6e-05 | 0.319295 ± 0.00011 | 0.299233 ± 0.0001 | 0.03413 |
| 100 | 0.220431 ± 7.9e-05 | 0.264368 ± 8.5e-05 | 0.242399 ± 8.1e-05 | 0.03354 |
| 110 | 0.198399 ± 5e-05 | 0.243610 ± 5.6e-05 | 0.221005 ± 5.3e-05 | 0.02408 |
| 120 | 0.179801 ± 5.1e-05 | 0.226065 ± 6.6e-05 | 0.202933 ± 5.8e-05 | 0.02875 |
| 130 | 0.163846 ± 3.6e-05 | 0.210986 ± 4.4e-05 | 0.187416 ± 4e-05 | 0.02138 |
| 140 | 0.150021 ± 2.5e-05 | 0.197890 ± 3.6e-05 | 0.173956 ± 3.1e-05 | 0.01755 |
| 150 | 0.137924 ± 6.7e-05 | 0.186430 ± 7.7e-05 | 0.162177 ± 7.2e-05 | 0.04424 |
| 160 | 0.127222 ± 3.2e-05 | 0.176283 ± 4.2e-05 | 0.151753 ± 3.7e-05 | 0.02436 |
| 170 | 0.117713 ± 6e-05 | 0.167279 ± 7.4e-05 | 0.142496 ± 6.7e-05 | 0.04715 |
| 180 | 0.109196 ± 3.1e-05 | 0.159237 ± 3.6e-05 | 0.134217 ± 3.4e-05 | 0.02515 |
| 190 | 0.101491 ± 1.5e-05 | 0.151963 ± 1.7e-05 | 0.126727 ± 1.6e-05 | 0.01273 |
| 200 | 0.094515 ± 2.2e-05 | 0.145400 ± 2.9e-05 | 0.119957 ± 2.5e-05 | 0.02116 |

## 1 kHz end checks

| Sweep | T_lower drift | T_upper drift | T_mean drift |
|---:|---:|---:|---:|
| 1 | +0.000770914 (+0.07742%) | +0.000700129 (+0.07046%) | +0.000735522 (+0.07394%) |
| 2 | -0.000554761 (-0.05564%) | -0.000499608 (-0.05022%) | -0.000527185 (-0.05293%) |
| 3 | +0.000716148 (+0.0719%) | +0.000704065 (+0.07084%) | +0.000710107 (+0.07137%) |
| 4 | -7.59101e-05 (-0.007618%) | -6.88424e-06 (-0.0006924%) | -4.13971e-05 (-0.004159%) |
| 5 | +0.000193359 (+0.01941%) | +0.000224238 (+0.02255%) | +0.000208798 (+0.02098%) |

## Anomalies and suitability

- Abnormal disagreement rule: a sweep point is flagged if it differs from its five-sweep T_mean by more than 3× the prior dataset maximum relative SD (0.1977%).
- No sweep point exceeded that threshold.
- Suitability: **suitable for IA analysis** based on completed sweeps, end checks, clipping checks, and recorded sideband SNR. Interpret the reported SDs as between-sweep repeatability; sideband SNR is estimated locally from FFT bins.

## Outputs

- `processed/formal_am_raw.csv` — one row per formal acquisition.
- `processed/formal_am_summary.csv` — per-frequency means and sample SDs across five sweeps.
- `processed/formal_am_end_checks.csv` — within-sweep initial-to-end drift.
- Native waveform archives are retained privately and are not included in this repository.
- `figures/` — requested plots.
