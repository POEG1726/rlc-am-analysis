# Drive-amplitude characterization summary

Date: 2026-09-26. Circuit: 4.7 nF polypropylene film capacitor, nominal 22 µH molded power inductor, and 5.1 Ω external resistor. The generator supplied a pure sine on CH2 with AM OFF. The sweep used 19 programmed amplitudes from 0.1 to 5.0 Vpp in ascending order and then descending order. Scope vertical ranges were selected automatically at each amplitude. At each setting, the phase-zero resonance, peak transfer, -3 dB bandwidth/Q, complex response, resonance Vin/Vout/current, and harmonics were recorded.

## Answer to the original issue

The original problem was that the measured phase-zero resonance moved with drive level (about 497.15 kHz at 1.2 Vpp versus 493.89 kHz at 5 Vpp), and peak transfer H also changed. That makes the circuit response amplitude-dependent under the original measurement conditions.

A real, repeatable amplitude dependence remains after correcting the scope range effect: measured f0 falls from about 498.40 kHz at 0.1 Vpp to 493.89 kHz at 5 Vpp; corrected Hmax falls from 0.723 to 0.682. The strongest evidence points to a circuit nonlinearity, with the molded power inductor the leading suspect, but the component source is not proven by this sweep alone.

The investigation also found a measurement-chain contribution: at very small signals, CH2 gain depended on the selected vertical range. For the 2 mV/div range, the scope read about 7.9% less amplitude than the 200 mV/div reference mapping, while CH1 range variation was around 0.5% or less. This created an artificial low-drive H increase in the raw data. Applying empirical relative range corrections brought Hmax at 0.1–0.2 Vpp into agreement (about 0.723 throughout). The correction removes the spurious low-drive H step; it does not remove the measured f0 trend.

## Findings

- **Stable low-drive region:** 0.1–0.2 Vpp is the flattest directly measured region. Across those settings, f0 spans about 62 Hz (0.012%), corrected Hmax about 0.0007 (0.10%), bandwidth about 116 Hz (0.23%), and Q about 0.024 (0.24%). The 0.1-to-0.2 Vpp mean resonance current rises from 0.87 to 1.74 mA peak. The complex RLC fit is not sufficiently accurate here (relative complex residual about 2.5%), so fitted Rs/L/C should not be used to claim those parameters have reached a plateau.
- **Onset:** small monotonic changes are already measurable by 0.3–0.6 Vpp: f0 has moved about 0.17–0.48 kHz below its 0.1 Vpp value. The departure becomes clear above roughly 1 Vpp. At 5 Vpp, relative to 0.1 Vpp, f0 is lower by 4.51 kHz (0.91%), Hmax by 5.77%, bandwidth is wider by 4.65%, and Q is lower by 5.31%.
- **Fitted parameters and resonance shift:** the complex fit meets the report's 2% residual criterion from 1.0 Vpp upward; lower-drive fits are marked unreliable. Over 1–5 Vpp, fitted Rs rises 3.27% (7.26 to 7.49 Ω), L_eff falls 1.59% (22.98 to 22.61 µH), and C_eff rises 3.10% (4.455 to 4.593 nF). Among the fitted parameters, C_eff has the largest fractional shift, with Rs close behind. The downward f0 shift mathematically follows the rise in the fitted L_eff·C_eff product. These are effective parameters from a simplified series RLC model; a fitted C_eff shift is not evidence that the polypropylene capacitor itself changed. The model can absorb frequency-dependent loss, parasitics, and inductor behavior into these fit values.
- **Harmonics:** harmonic levels generally rise with drive. Between 0.2 and 5 Vpp, CH1 H2/H3 rise by about 19/31 dB, and CH2 H2/H3 by about 28/18 dB. CH1 already contains notable harmonics, so output harmonics cannot be attributed to the inductor without separately characterizing the generator and circuit. At the highest setting the measured CH1 H2/H3 are about -44/-39 dBc; CH2 values are about -52/-57 dBc.
- **Heating or hysteresis:** ascending and descending sweeps agree closely. Maximum paired differences are 52 Hz in f0, 0.00184 in Hmax, 123 Hz in bandwidth, and 0.023 in Q. This gives no clear evidence of thermal drift or hysteresis over this sweep and timing.
- **Likely source:** the inductor remains the leading physical suspect because amplitude-dependent inductance/loss in a molded power inductor is plausible, while a polypropylene film capacitor is expected to be comparatively linear. However, this experiment did not independently measure the inductor or capacitor, and the simplified fit does not isolate component behavior. The evidence supports “likely inductor or other circuit nonlinearity,” not a proven attribution.
- **IA drive recommendation:** use **0.2 Vpp programmed** as the formal drive starting point. It lies inside the measured 0.1–0.2 Vpp flat region and gives about 24.6 mVpp at measured CH1 Vin, 8.9 mV peak at CH2 Vout, and 1.74 mA peak through 5.1 Ω. This is a reasonable SNR level for the existing long-record scope capture. Before the full AM-frequency dataset, make one pilot modulation capture and confirm sideband SNR at the chosen modulation depth; if the sidebands are not cleanly resolved, test 0.3 Vpp and recheck f0/H stability. No AM modulation-frequency dataset was collected in this characterization.

## Scope calibration caveat

The scale correction was derived by measuring the same sine while stepping adjacent vertical ranges, discarding the first capture after each change, and averaging three captures per range. Factors map each range to the 200 mV/div reference; they correct relative range gain only and are not an absolute scope calibration. The corrected actual Vin/Vout/H/current values are therefore better matched across auto-ranged settings than the raw values, but retain scope calibration uncertainty. Frequency, phase, bandwidth, Q, and dBc harmonic ratios do not depend on this amplitude correction.

## Deliverables

- `drive_amplitude_sweep.csv`: one row per amplitude and sweep direction; raw resonance and harmonic columns are retained with `_uncorrected` suffixes where corrected values replace them. Includes relative range correction factors and refitted effective parameters.
- `drive_amplitude_frequency_scan.csv`: all coarse/fine/BW frequency samples, raw and range-corrected complex transfer data.
- `drive_amplitude_harmonics.csv`: H2/H3 values with corrected absolute line amplitudes.
- `vertical_scale_corrections.json` and `scope_vertical_range_calibration.csv`: empirical range correction method and calibration captures.
- `drive_amplitude_log.txt`: instrument settings, sweep progress, and final restored state.

Plots: `f0_vs_programmed_vpp.png`, `f0_vs_actual_vin.png`, `f0_vs_resonance_current.png`, `hmax_vs_programmed_vpp.png`, `hmax_vs_current.png`, `fitted_parameters_vs_amplitude.png`, `fitted_parameters_vs_current.png`, `bandwidth_q_vs_amplitude.png`, `harmonics_vs_amplitude.png`, and `ascending_descending_comparison.png`.

The generator was restored to CH2 sine, 5 Vpp at 493.9 kHz with AM OFF and output ON; the scope's original 50 ms/div acquisition state and channel vertical ranges were restored after measurement.
