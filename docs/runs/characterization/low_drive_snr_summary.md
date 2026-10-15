# Low-drive AM sideband SNR and repeatability

Date: 2026-09-26. Operating point: generator CH2 sine carrier at 0.2 Vpp, 498.3374 kHz programmed (approximately the low-drive phase-zero resonance); internal sine AM at 50% depth. CH1 measured the actual RLC input and CH2 the resistor output. The sweep covered 1, 2, 3, 5, 8, 10, 12, 15, 18, 20, 22, 25, 30, 35, and 40 kHz, with three repeated captures at each frequency.

## Acquisition and analysis

- Each capture used 5,000,000 native waveform points per channel, WAV interval 1 (no decimation), 50 MSa/s, 0.100 s record, and 10 Hz FFT-bin spacing. Scope ranges stayed fixed for all repeats and frequencies: CH1 20 mV/div and CH2 10 mV/div.
- For each channel, the Hann-windowed FFT peak was searched within ±100 Hz of each expected LSB/carrier/USB. H was calculated from corrected CH2/CH1 line amplitudes; T_lower=H_L/H_C, T_upper=H_U/H_C, and T_mean=(T_lower+T_upper)/2. A fixed per-channel scale correction affects absolute H but cancels in normalized T and same-channel SNR.
- Local noise floor is the median FFT-bin amplitude in the two sideband-adjacent windows 100–450 Hz from the target line. Bins more than 12 dB above the preliminary median were treated as narrow spurs and excluded with a ±20 Hz guard. SNR is the peak-bin amplitude divided by that local median-bin floor, in dB. This is a local spectral-bin SNR estimate, not an integrated broadband noise measurement.

## Results

- Worst measured sideband SNR was **84.5 dB**, for **CH1 USB at 1 kHz**. The minimum across each of the four channel/sideband combinations was: CH1 LSB 86.3 dB; CH1 USB 84.5 dB; CH2 LSB 87.3 dB; CH2 USB 85.4 dB.
- Repeatability was strong within this sweep: maximum relative sample SD was **{max_l:.3f}% for T_lower at {fm_l/1000:g} kHz** and **{max_u:.3f}% for T_upper at {fm_u/1000:g} kHz** (n=3). The largest absolute SD was {max(s.T_lower_sd.max(),s.T_upper_sd.max()):.6f}.
- No capture clipped. The largest observed code magnitude used {100*maxocc:.1f}% of the signed ADC code range. For every fm, CH1 and CH2 selected the same FFT bin for each of LSB/carrier/USB, and all six line-bin indices were stable across the three repeats.
- T decreases at higher fm because the sidebands move farther from resonance: mean T_mean is 0.9971 at 1 kHz and 0.5351 at 40 kHz. This is the measured circuit response, while the sideband SNR remains high across the sweep.

| fm (kHz) | T_lower mean ± SD | T_upper mean ± SD | T_mean | Worst SNR (dB) | max relative SD of T (%) |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.997871 ± 0.000106 | 0.996346 ± 0.000199 | 0.997109 | 84.5 | 0.020 |
| 2 | 0.996186 ± 0.000140 | 0.993556 ± 0.000162 | 0.994871 | 85.4 | 0.016 |
| 3 | 0.993099 ± 0.000170 | 0.989407 ± 0.000161 | 0.991253 | 85.7 | 0.017 |
| 5 | 0.982399 ± 0.000340 | 0.976909 ± 0.000326 | 0.979654 | 87.9 | 0.035 |
| 8 | 0.956403 ± 0.000473 | 0.949379 ± 0.000484 | 0.952891 | 88.1 | 0.051 |
| 10 | 0.933644 ± 0.000305 | 0.926589 ± 0.000335 | 0.930116 | 88.5 | 0.036 |
| 12 | 0.907198 ± 0.000338 | 0.900763 ± 0.000324 | 0.903981 | 89.2 | 0.037 |
| 15 | 0.863341 ± 0.000072 | 0.858857 ± 0.000085 | 0.861099 | 87.0 | 0.010 |
| 18 | 0.816393 ± 0.000301 | 0.815026 ± 0.000272 | 0.815709 | 89.3 | 0.037 |
| 20 | 0.784757 ± 0.000137 | 0.785881 ± 0.000165 | 0.785319 | 87.6 | 0.021 |
| 22 | 0.753824 ± 0.000180 | 0.757440 ± 0.000187 | 0.755632 | 88.2 | 0.025 |
| 25 | 0.708690 ± 0.000150 | 0.716192 ± 0.000198 | 0.712441 | 89.4 | 0.028 |
| 30 | 0.638733 ± 0.000172 | 0.652474 ± 0.000183 | 0.645604 | 88.7 | 0.028 |
| 35 | 0.576833 ± 0.000079 | 0.596065 ± 0.000074 | 0.586449 | 88.5 | 0.014 |
| 40 | 0.523156 ± 0.000133 | 0.547061 ± 0.000134 | 0.535108 | 88.6 | 0.025 |

## Adequacy for the formal IA

Using a practical screening criterion of at least 20 dB local peak-bin SNR on both channels for both sidebands and no more than 1% relative SD in T across the three captures, **0.2 Vpp is adequate across 1–40 kHz**. The weakest sideband clears the SNR criterion by more than 64 dB, and the worst within-sweep T relative SD is about 0.051%. No higher drive amplitude needs to be tested for SNR on this evidence.

One earlier low-drive 1 kHz single-capture check gave T_mean=0.99968, versus 0.99711 for the present three-capture mean (difference 0.00257, or 0.26%). Both are close to unity, and the present within-sweep repeats are much tighter. Treat this as a small between-setup/run offset; keep the carrier, scope scales, and measurement procedure fixed for the formal dataset and retain repeated captures.

## Deliverables

- `low_drive_snr_sweep.csv`: aggregated one-row-per-modulation-frequency results, means/SDs, minimum SNR, clipping and bin-consistency flags.
- `low_drive_snr_repeats.csv`: all 45 per-capture line amplitudes, H/T, noise floors, SNRs, bins, and ADC extrema.
- `low_drive_snr.png`: sideband SNR by channel and sideband (mean ± 1 SD).
- `low_drive_sideband_repeatability.png`: T means with ±1 SD and relative repeatability versus fm.
- `low_drive_snr_sweep.log`: instrument settings, capture progress, and final restoration readback.

The generator was restored to CH2 sine, 5 Vpp at 493.9 kHz, AM OFF, HiZ, output ON. The scope was restored to its prior 50 ms/div, 50 MSa/s, 25 Mpoint state and 200/100 mV/div ranges. No formal AM modulation-frequency dataset was collected.
