# Low-drive 1 kHz AM verification

Date: 2026-09-26. Repeated the original single-point AM diagnostic at the low-drive operating point identified by the amplitude sweep:

- Generator CH2 carrier: sine, 0.2 Vpp programmed, 498.3374 kHz programmed (the mean low-drive phase-zero setting; measured carrier peak 498.3434 kHz).
- AM: internal sine, 1.000 kHz, 50% depth; output load HiZ; output ON.
- Scope: CH1 measured RLC input, CH2 resistor output; 1 MΩ inputs; 50 MSa/s; paired 5,000,000-point native captures (0.100 s, 10 Hz FFT spacing). Fixed ranges were 20 mV/div on CH1 and 10 mV/div on CH2 to retain headroom at the AM envelope peaks.

## Measured result

| Component | Measured peak frequency | Corrected Vin (Vpk) | Corrected Vout (Vpk) | Corrected H |
|---|---:|---:|---:|---:|
| Lower sideband | 497,343.43 Hz | 0.002876 | 0.002073 | 0.720586 |
| Carrier | 498,343.44 Hz | 0.011490 | 0.008275 | 0.720204 |
| Upper sideband | 499,343.45 Hz | 0.002878 | 0.002070 | 0.719362 |

- **T_lower = 1.00053**
- **T_upper = 0.99883**
- **T_mean = 0.99968**

The three peaks are about 6 Hz above programmed values, consistent with the small generator/scope frequency-axis offset seen in previous captures. Both channels selected the same FFT bin for each component. No ADC clipping was present (CH1 code range -6784 to 7632; CH2 -10480 to 10400, versus signed limits ±32768).

The range correction changes absolute H slightly; the within-channel normalized T values are insensitive to a fixed range gain because it cancels between sideband and carrier. A pure-sine low-drive interpolation at these same three frequencies gives T_lower ≈ 0.99934, T_upper ≈ 0.99878, and T_mean ≈ 0.99906. The AM measurement agrees closely with that linear small-signal response and with the expected approximately 0.9992 value for Q≈9.8 at 1 kHz offset.

## Interpretation

This low-drive, resonance-tuned AM point does **not** reproduce the earlier 0.943 normalized mean. It gives approximately 1.000. Together with the measured amplitude-dependent resonance and transfer curve, this strongly supports high-drive amplitude-dependent circuit behavior as the main reason the earlier 1 kHz sidebands were attenuated. It does not identify the exact component mechanism; the molded power inductor remains a leading suspect, not a proven source. This was one verification capture, not an AM modulation-frequency sweep.

## Files

- `am_1khz_low_drive.csv`: line-by-line corrected and raw FFT amplitudes, gains, peak frequencies, and bins.
- `am_1khz_low_drive_summary.csv`: one-row conditions and normalized result.
- `am_1khz_low_drive.npz`: paired native waveform vectors and metadata.
- `am_1khz_low_drive_spectrum.png`: spectrum around the AM triplet.
- `am_1khz_low_drive.log`: instrument IDs, readbacks, acquisition and restored state.

After acquisition the instruments were restored to CH2 sine 5 Vpp at 493.9 kHz, AM OFF, output ON/HiZ, and the scope's prior 50 ms/div timebase and CH1/CH2 ranges.
