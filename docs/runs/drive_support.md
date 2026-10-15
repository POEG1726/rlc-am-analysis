# Fresh drive-dependence supporting measurement

- Run completed 2026-10-01 UTC. This summary uses only the 18 acquisitions collected in this run; no earlier measurement-result files were loaded.
- AM was OFF throughout. Existing pure-sine workflow used a coarse response scan, phase-zero refinement, Hmax and -3 dB crossings, complex transfer fitting, and a resonance waveform for Vin/Vout/current and harmonic checks.
- The scope was set to 50 MSa/s, 100,000 waveform points, 1 MΩ inputs, DC coupling, and 1× probes. Vertical scales were auto-selected from live waveforms for every repeat.
- The acquisition wait was increased to 500 ms after live diagnostics found the scope could return an empty waveform block after a 120 ms wait. A separate diagnostic run verified 10/10 full-length reads at 500 ms; those diagnostics were not included in these results.

## Results

- **Acquisitions:** 18/18 (0.2, 1, 2, 3, 4, 5 Vpp × 3 repeats). **All 18 passed** validity checks.
- **Clipping:** none in resonance captures or any of the 918 frequency-scan captures.
- Phase-zero resonance moved from **498.004 kHz at 0.2 Vpp** to **493.519 kHz at 5 Vpp** (Δf₀ = **-4.485 kHz**, -0.901%).
- The first measured step already shows a shift: 0.2 to 1 Vpp changes mean f₀ by **-860 Hz**. Since amplitudes between those settings were not tested in this supporting run, the onset is bounded to **above 0.2 and at or below 1 Vpp**, not pinpointed.
- Hmax changed by **-0.06481** (-8.65%) from 0.2 to 5 Vpp. Mean bandwidth changed from 50.69 to 53.04 kHz and Q from 9.825 to 9.304.
- Resonance current, calculated from Vout,pk / 5.1 Ω, rose from 1.800 mA peak at 0.2 Vpp to 43.684 mA peak at 5 Vpp.

| Drive (Vpp) | f₀ mean (kHz) | f₀ sample SD (Hz) | Hmax mean | BW mean (kHz) | Q mean | Vin at f₀ (Vpp) | Vout at f₀ (Vpk) | Ires (mA pk) | CH1/CH2 scale (V/div) | Valid repeats |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|
| 0.2 | 498.00383 | 1.48 | 0.749476 | 50.686 | 9.8253 | 0.0245124 | 0.00918225 | 1.8004 | 0.01/0.005 | 3/3 |
| 1 | 497.14407 | 0.58 | 0.717637 | 51.319 | 9.6873 | 0.123844 | 0.0444275 | 8.7113 | 0.05/0.02 | 3/3 |
| 2 | 496.14626 | 0.53 | 0.704137 | 51.887 | 9.5620 | 0.251062 | 0.0883831 | 17.3300 | 0.1/0.05 | 3/3 |
| 3 | 495.25109 | 0.62 | 0.696832 | 52.289 | 9.4715 | 0.379982 | 0.132385 | 25.9578 | 0.1/0.1 | 3/3 |
| 4 | 494.27093 | 0.52 | 0.690798 | 52.765 | 9.3674 | 0.517686 | 0.17876 | 35.0510 | 0.2/0.1 | 3/3 |
| 5 | 493.51897 | 1.52 | 0.684669 | 53.043 | 9.3041 | 0.65083 | 0.222786 | 43.6836 | 0.2/0.1 | 3/3 |

## Interpretation

- **Stable low-drive point:** 0.2 Vpp was repeatable in these three captures (f₀ sample SD 1.48 Hz; Hmax sample SD 0.000601) and gave the highest measured Hmax. It is the recommended formal-IA drive among the tested settings, with the chosen 10 mV/div and 5 mV/div scales providing usable resonance waveforms without clipping.
- **Amplitude dependence:** clearly present by 1 Vpp and grows through 5 Vpp: f₀ and Hmax decrease while bandwidth increases and Q decreases. The data do not locate a sharper threshold between 0.2 and 1 Vpp.
- **Complex-fit parameters:** all three fits at 0.2 Vpp failed the workflow reliability criterion (complex-H RMSE about 0.034); do not use them to claim low-drive L/C/R values. For the reliable 1–5 Vpp fits, mean Rs changes 0.223 Ω (+3.09%), L_eff changes -0.924 µH (-4.00%), and C_eff changes 0.257 nF (+5.79%). C_eff has the largest relative shift in that fitted range (5.79%). These are correlated equivalent-circuit fit parameters; the fitted C_eff shift should not be read as proof the polypropylene capacitor physically changed.
- **Harmonics:** mean CH2 H2 changes from -76.2 dBc at 0.2 Vpp to -51.0 dBc at 5 Vpp; H3 changes from -77.3 to -57.4 dBc. These rises are consistent with more nonlinear distortion at higher drive; the very low-drive harmonic values may be influenced by the measurement floor.
- **Repeatability / heating / hysteresis:** the within-amplitude f₀ range across three repeats was at most 2.99 Hz, and no clipping or failed captures occurred. Because this run used increasing amplitudes only and did not include a descending sweep or temperature measurement, it cannot establish absence of thermal drift or hysteresis.
- **Likely component source:** the current-dependent resonance and harmonic changes make the molded power inductor a plausible contributor, especially given its magnetic core. This experiment measured the assembled circuit and does not isolate the inductor from the generator, resistor, wiring, or fit-model limitations; the source is therefore **not proven**. The fit trend is not a unique component-level diagnosis.

## Deliverables

- `drive_dependence_support.csv` — 18 individual fresh resonance results, including measured Vin/Vout, calculated peak/RMS current, per-repeat scope scales, clipping/validity, bandwidth/Q, fit parameters and harmonics.
- `drive_dependence_frequency_scans.csv` — 918 per-frequency records from coarse/refined scans and direct resonance captures.
- `drive_dependence_harmonics.csv` — 2nd/3rd harmonic checks at each resonance.
- `f0_vs_drive.png`, `Hmax_vs_drive.png`, `f0_vs_current.png` — requested plots with individual repeats and means/sample SD where applicable.
- `drive_dependence_support.log` — instrument identities/readbacks, scope settings and per-acquisition progress.

No AM-frequency dataset was collected.
