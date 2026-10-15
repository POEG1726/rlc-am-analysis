# Fresh drive-dependence data with 12 Vpp extreme test

- Run timestamp: 2026-10-01 (UTC date; time withheld).
- Run folder: `drive_dependence_extreme_12Vpp_20260930_run02`. Recovery attempt: `resume_4V_repeat2_after_IO_fault_01`; the incomplete 4 Vpp repeat 2 scan remains preserved in raw records and is excluded from the final response overlay.
- AM OFF, CH2 sine, 50 MSa/s, 100,000 points/channel, WAV interval 1, DC coupling, 1× probes, 1 MΩ scope inputs; live auto-ranging was performed from this run’s fresh waveforms.
- Run status: **COMPLETE**; valid completed repeats: **21/21**.
- Clipping in saved result captures: **no**.
- Safety/validity stop reason: **none**.

## VISA communication interruption and recovery

- The original 4 Vpp repeat 2 coarse scan stopped at 545 kHz with `VI_ERROR_IO (-1073807298)`; the failed transfer was not treated as a measurement.
- A link check then passed 30 `*IDN?` queries, 20 full native waveform channel transfers, and 3 fresh 4 Vpp / 545 kHz captures. A separate first-capture attempt in `run03` later produced the same generic I/O error at 0.2 Vpp and was stopped with the generator off.
- The run resumed at 4 Vpp repeat 2, restarting that repeat from its initial range check. All remaining measurements then completed with valid records.
- This pattern supports an intermittent VISA/acquisition-path fault. It does not establish whether the physical USB connection caused it.

| Drive (Vpp) | f₀ mean (kHz) | f₀ SD (Hz) | Hmax | BW (kHz) | Q | Vin at f₀ (Vpp) | Vout at f₀ (Vpk) | Ires (mA pk) | reliable fits | valid repeats |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 497.98073 | 1.23 | 0.749698 | 50.659 | 9.8301 | 0.0244914 | 0.00917911 | 1.7998 | 0/3 | 3/3 |
| 1 | 497.12250 | 1.62 | 0.717591 | 51.337 | 9.6836 | 0.123815 | 0.0444052 | 8.7069 | 3/3 | 3/3 |
| 2 | 496.12259 | 1.92 | 0.704513 | 51.909 | 9.5575 | 0.250895 | 0.0883705 | 17.3275 | 3/3 | 3/3 |
| 3 | 495.23354 | 0.706 | 0.696552 | 52.305 | 9.4682 | 0.380161 | 0.132386 | 25.9580 | 3/3 | 3/3 |
| 4 | 494.25432 | 11 | 0.691519 | 52.638 | 9.3897 | 0.51728 | 0.178828 | 35.0643 | 3/3 | 3/3 |
| 5 | 493.51098 | 2.25 | 0.685952 | 52.883 | 9.3322 | 0.650091 | 0.222939 | 43.7135 | 3/3 | 3/3 |
| 12 **— EXTREME DRIVE TEST** | 489.23397 | 2.64 | 0.650839 | 54.490 | 8.9784 | 1.62301 | 0.528193 | 103.5673 | 3/3 | 3/3 |

- 0.2 to 12 Vpp f₀ difference: **-8.7468 kHz**.
- 0.2 to 12 Vpp Hmax difference: **-0.09886**.
- 12 Vpp mean resonance current: **103.567 mA peak**.

## Files

- `drive_dependence_results.csv` — per-repeat phase-zero resonance, Hmax, bandwidth/Q, Vin/Vout/current, fit parameters, harmonics, ranges and clipping/validity fields.
- `frequency_response_raw.csv` — every fresh complex H frequency point, saved incrementally; includes the complete resonance response curves.
- `acquisition_audit.csv` — per-capture point count, code extrema, clipping, generator/scope status, and validity.
- `harmonics.csv` — individual repeat H2/H3 measurements on CH1 and CH2.
- `plots/` — every requested plot, with 12 Vpp in red and explicitly labelled.

12 Vpp is an extreme diagnostic point and is not recommended as a normal operating condition.
