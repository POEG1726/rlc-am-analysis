# RLC characterization summary

## Result

The physical circuit is resonant near **493.887 kHz** at the tested 5 Vpp generator setting. Use **493.9 kHz programmed** as the starting carrier under this same drive condition. AM remains OFF. The 500 kHz setting was about 6.1 kHz above the phase-zero resonance.

The broad 450–550 kHz pure-sine sweep found a broad response, not an unusually narrow passband:

- Maximum measured |Vout/Vin|: **0.6827** (the peak is very flat; the sampled maximum is at 494.211 kHz).
- Phase-zero resonance from the 100 Hz refinement: **493887 Hz programmed** (linear phase fit RMS 0.002°).
- Measured -3 dB frequencies: **467.953 kHz** and **521.047 kHz**.
- Measured bandwidth: **53.09 kHz**, Q = **9.30**.
- Complex impedance fit: total series resistance **7.496 Ω**; effective inductance **22.609 µH**; effective capacitance **4.593 nF**; model resonance **493.903 kHz**; model bandwidth **52.76 kHz** and Q **9.36**.
- Subtracting the nominal 5.1 Ω output resistor gives about **2.396 Ω** of other effective series loss. This includes inductor AC loss and PCB/capacitor losses; it cannot be assigned to the inductor alone. The separate 0.21 Ω inductor DCR is not a suitable AC-loss estimate.

The effective L/C values are fitted, not direct component measurements. The model assumes the nominal 5.1 Ω resistor and approximately constant series loss over the fit span. The fit's complex-H RMS residual is 0.0030; phase-derived delay is -13.0 ps. Values should therefore be reported at sensible precision rather than as exact component values.

## Diagnosis of the 1 kHz AM result

I repeated the specified 500 kHz, 1 kHz, 50% depth, 5 Vpp CH2 AM measurement with a 5 Mpoint native capture. It gave **T_lower=0.9538, T_upper=0.9308, T_mean=0.9423**, reproducing the reported 0.954/0.931/0.943 result. Both channels selected the same FFT bin for each line. The measured spectral peaks were approximately 5.6 Hz above their programmed values, consistent with a small generator/scope clock offset.

After retuning the carrier to 493.9 kHz, a separate AM diagnostic still gave T_lower=0.9686, T_upper=0.9575, T_mean=0.9631; the pure-sine curve there predicts a normalized mean close to 1.000. Retuning improves the result but does not remove the sideband deficit. Pure-sine measurements at the same frequencies predict H values of about 0.6700, 0.6647, and 0.6591; normalized to the 500 kHz pure-sine carrier, their mean is **0.9997**, not 0.943. At resonance, the pure-sine response is similarly flat over ±1 kHz. The carrier H from AM agrees with the pure-sine 500 kHz H within about 1%, while the AM sideband H values are lower by several percent.

A second physical effect was measured: reducing generator amplitude from 5 Vpp to 1.2 Vpp changed CH1 from roughly 0.33 Vpk to 0.075 Vpk, changed H near 500 kHz from about 0.665 to 0.709, and shifted the phase-zero frequency from about 493.89 kHz to 497.15 kHz. That is direct evidence that this operating circuit is amplitude-dependent. It is consistent with nonlinear or drive-dependent RLC/parasitic behavior, but these measurements do not identify which component causes it. Under 50% AM, the carrier and sidebands therefore do not behave as independent small-signal sine tests.

The discrepancy is **not explained by a narrow measured bandwidth, small carrier detuning, FFT-bin choice, transfer scaling, clipping, or the 0.21 Ω DCR estimate**. The most likely physical explanation is drive-level-dependent response under AM; the exact mechanism remains unlocalized. The pure-sine sweep and AM data both use CH1 as the measured input reference, so generator programmed amplitude is not used as Vin.

## Pipeline and acquisition checks

- Generator queried and controlled on **CH2 only**; CH2 output is ON, HiZ, sine, zero offset. AM was OFF during the sine sweep and is OFF at completion.
- Scope CH1/CH2 are DC coupled, 1× probes, 1 MΩ inputs (not 50 Ω); bandwidth limit is 20 MHz. Vertical scales were 200 mV/div and 100 mV/div.
- Scope waveform descriptor reported dt = 20 ns (50 MSa/s); transfer interval was 1 (no decimation). Sweep records used 100,000 points; AM diagnostics used 5,000,000 points (0.100 s and 10 Hz FFT spacing).
- No clipping: AM-capture voltage ranges were approximately CH1 [-0.527, 0.525] V and CH2 [-0.330, 0.338] V, within the configured ranges.
- Existing Hann/rFFT amplitude extraction and local peak search reproduced the supplied AM values. Its peak-bin selection is valid for these high-SNR components. Keep the actual FFT bins and channel-specific peaks in the CSV for audit; a complex sine fit is preferable for future pure-sine impedance measurements.
- Three repeated near-resonance pure-sine captures agreed in H within 0.022% standard deviation. The fresh 500 kHz AM result and the supplied prior result agree within 0.001 in each normalized sideband response; the 493.9 kHz AM point was acquired once.

## Readiness

The instrument control, raw waveform transfer, CH1-referenced FFT ratios, and resonance scan are validated. The setup is **not yet validated as a linear, amplitude-independent RLC experiment**. Before formal IA data collection, either hold the 5 Vpp/50% AM operating point fixed and frame the result as an empirical driven response, or lower/choose the drive level and confirm that the resonance and transfer curve stop shifting with amplitude. Do not interpret the 0.943 result as the linear small-signal response of a fixed passive RLC without that control.

No large modulation-frequency dataset was collected. Final instrument state: generator CH2 5 Vpp sine at 493.9 kHz, output ON/HiZ, AM OFF; scope returned to 50 ms/div, 50 MSa/s, 25 Mpoint acquisition with CH1/CH2 at 1 MΩ.

## Files

- `resonance_scan.csv`: coarse and fine sweep values, repeated points, complex H and derived impedance.
- `resonance_refinement.csv`: 100 Hz resonance refinement.
- `characterization_parameters.json`: measured and fitted parameter summary.
- `resonance_curve.png`, `impedance_real_imag.png`, `drive_level_comparison.png`: response, impedance, and drive-level comparison plots.
- `am_1khz_diagnostic_500khz.csv` and `.npz`, plus `am_1khz_diagnostic.csv` and `.npz`: 500 kHz and resonance AM line amplitudes with compressed paired waveform vectors.
- `amplitude_comparison.csv`, `amplitude_linearity_check.csv`: low-drive response and matching sideband-level pure-sine checks.
- `characterization_log.txt`, `am_1khz_diagnostic_500khz.log`: instrument IDs, settings, scan values, and diagnostic acquisition details.
- `characterize_rlc.py`, `refine_resonance.py`, `diagnose_am_1khz_500khz.py`, `compare_drive_levels.py`: acquisition and analysis scripts.


## Low-drive 1 kHz AM follow-up (2026-09-26)

A single repeat of the original 1 kHz, 50% AM test was acquired at 0.2 Vpp with the carrier tuned to the low-drive phase-zero setting (498.3374 kHz programmed; measured peak 498.3434 kHz). The normalized sideband results were T_lower=1.00053, T_upper=0.99883, T_mean=0.99968. A low-drive pure-sine interpolation predicts T_mean≈0.99906 at the same three frequencies, in close agreement. The 5 Mpoint paired capture showed no ADC clipping, both channels selected the same FFT bins, and line frequencies were consistently about 6 Hz above their programmed values.

This low-drive result does not reproduce the earlier T_mean≈0.943. Along with the drive-amplitude sweep, it supports amplitude-dependent circuit behavior at the earlier high drive as the main explanation for the AM sideband deficit. It does not prove which component causes the nonlinearity; the power inductor remains a leading suspect. This was one AM verification point; no modulation-frequency dataset was collected. See `am_1khz_low_drive_summary.md` for details and files.
