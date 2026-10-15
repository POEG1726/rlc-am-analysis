
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import pyvisa
import numpy as np
import matplotlib.pyplot as plt
import struct
import re

SCOPE_ADDR = _bench.scope_address
GEN_ADDR   = _bench.generator_address

# First test:
# 5 Mpts @ 50 MSa/s = 0.1 s -> 10 Hz FFT bin spacing
POINTS_TO_READ = 5_000_000


rm = pyvisa.ResourceManager()

scope = rm.open_resource(SCOPE_ADDR)
gen = rm.open_resource(GEN_ADDR)

scope.timeout = 120_000
gen.timeout = 5_000

# Large USB transfers
scope.chunk_size = 16 * 1024 * 1024


def query_binary_block(inst, command):
    inst.write(command)
    raw = inst.read_raw()

    pos = raw.find(b"#")
    if pos < 0:
        raise RuntimeError(
            f"No IEEE binary block found.\n"
            f"Command: {command}\n"
            f"Beginning of response: {raw[:100]!r}"
        )

    n_digits = int(chr(raw[pos + 1]))

    length_start = pos + 2
    length_end = length_start + n_digits

    data_length = int(raw[length_start:length_end].decode())

    data_start = length_end
    data_end = data_start + data_length

    block = raw[data_start:data_end]

    if len(block) != data_length:
        raise RuntimeError(
            f"Expected {data_length} binary bytes, "
            f"received {len(block)}."
        )

    return block


def read_channel(ch, points):
    scope.write(f":WAV:SOUR C{ch}")

    # SDS824X HD 12-bit acquisition -> WORD transfer
    scope.write(":WAV:WIDT WORD")
    scope.write(":WAV:BYT LSB")

    # IMPORTANT: no transfer decimation
    scope.write(":WAV:INT 1")

    scope.write(":WAV:STAR 0")
    scope.write(f":WAV:POIN {points}")

    pre = query_binary_block(scope, ":WAV:PRE?")

    vdiv = struct.unpack_from("<f", pre, 156)[0]
    voffset = struct.unpack_from("<f", pre, 160)[0]
    code_per_div = struct.unpack_from("<f", pre, 164)[0]
    dt = struct.unpack_from("<f", pre, 176)[0]

    data = query_binary_block(scope, ":WAV:DATA?")

    codes = np.frombuffer(data, dtype="<i2")

    voltage = (
        codes.astype(np.float64)
        * (vdiv / code_per_div)
        - voffset
    )

    return voltage, dt


def get_generator_frequencies():
    s = gen.query("C2:MDWV?").strip()

    print("\nGenerator modulation state:")
    print(s)

    # Internal modulation frequency
    m = re.search(
        r"SRC,INT,FRQ,([0-9.+\-Ee]+)HZ",
        s
    )

    # Carrier frequency
    c = re.search(
        r"CARR,.*?FRQ,([0-9.+\-Ee]+)HZ",
        s
    )

    if not m or not c:
        raise RuntimeError(
            "Could not parse carrier/modulation frequency "
            "from C2:MDWV?"
        )

    fm = float(m.group(1))
    fc = float(c.group(1))

    return fc, fm


def fft_spectrum(x, fs):
    # Remove any DC offset
    x = x - np.mean(x)

    n = len(x)

    # Hann window
    window = np.hanning(n)
    coherent_gain = np.mean(window)

    y = np.fft.rfft(x * window)

    freqs = np.fft.rfftfreq(n, d=1/fs)

    # Single-sided peak-amplitude spectrum
    amp = 2.0 * np.abs(y) / (n * coherent_gain)

    # DC is not doubled
    amp[0] /= 2

    return freqs, amp


def peak_near(freqs, amp, target, search_width=100):
    mask = (
        (freqs >= target - search_width)
        & (freqs <= target + search_width)
    )

    if not np.any(mask):
        raise RuntimeError(
            f"No FFT bins around {target} Hz"
        )

    local_freqs = freqs[mask]
    local_amp = amp[mask]

    i = np.argmax(local_amp)

    return local_freqs[i], local_amp[i]


print("Scope:")
print(scope.query("*IDN?").strip())

print("\nGenerator:")
print(gen.query("*IDN?").strip())

print("\nScope settings:")
print("Acquisition points:", scope.query(":ACQ:POIN?").strip())
print("Sample rate:", scope.query(":ACQ:SRAT?").strip())


fc, fm = get_generator_frequencies()

print(f"\nCarrier fc     = {fc:,.3f} Hz")
print(f"Modulation fm  = {fm:,.3f} Hz")

fl = fc - fm
fu = fc + fm

print(f"Expected LSB    = {fl:,.3f} Hz")
print(f"Expected carrier= {fc:,.3f} Hz")
print(f"Expected USB    = {fu:,.3f} Hz")


# Freeze one acquisition so CH1 and CH2 come from same record
scope.write(":TRIG:STOP")

print("\nReading CH1...")
vin, dt1 = read_channel(1, POINTS_TO_READ)

print("Reading CH2...")
vout, dt2 = read_channel(2, POINTS_TO_READ)

scope.write(":TRIG:RUN")


if not np.isclose(dt1, dt2):
    raise RuntimeError(
        f"CH1 dt = {dt1}, CH2 dt = {dt2}"
    )

fs = 1.0 / dt1

print("\nDownloaded waveform:")
print(f"Points          = {len(vin):,}")
print(f"dt              = {dt1:.12g} s")
print(f"Sample rate     = {fs/1e6:.6f} MSa/s")
print(f"Record duration = {len(vin)/fs:.6f} s")
print(f"FFT bin spacing = {fs/len(vin):.6f} Hz")


print("\nCalculating CH1 FFT...")
freq, amp_in = fft_spectrum(vin, fs)

print("Calculating CH2 FFT...")
_, amp_out = fft_spectrum(vout, fs)


# Extract the three AM components
fin_l, vin_l = peak_near(freq, amp_in, fl)
fin_c, vin_c = peak_near(freq, amp_in, fc)
fin_u, vin_u = peak_near(freq, amp_in, fu)

fout_l, vout_l = peak_near(freq, amp_out, fl)
fout_c, vout_c = peak_near(freq, amp_out, fc)
fout_u, vout_u = peak_near(freq, amp_out, fu)


HL = vout_l / vin_l
HC = vout_c / vin_c
HU = vout_u / vin_u

# Normalize sideband transmission against carrier transmission
T_lower = HL / HC
T_upper = HU / HC
T_mean = (T_lower + T_upper) / 2


print("\n================ FFT RESULTS ================")

print("\nCH1 — RLC input")
print(f"LSB     {fin_l:12.2f} Hz   {vin_l:.6f} Vpk")
print(f"Carrier {fin_c:12.2f} Hz   {vin_c:.6f} Vpk")
print(f"USB     {fin_u:12.2f} Hz   {vin_u:.6f} Vpk")

print("\nCH2 — resistor output")
print(f"LSB     {fout_l:12.2f} Hz   {vout_l:.6f} Vpk")
print(f"Carrier {fout_c:12.2f} Hz   {vout_c:.6f} Vpk")
print(f"USB     {fout_u:12.2f} Hz   {vout_u:.6f} Vpk")

print("\nTransfer magnitude")
print(f"H_L = {HL:.6f}")
print(f"H_C = {HC:.6f}")
print(f"H_U = {HU:.6f}")

print("\nNormalized sideband transmission")
print(f"T_lower = {T_lower:.6f}")
print(f"T_upper = {T_upper:.6f}")
print(f"T_mean  = {T_mean:.6f}")


# -------- Plot only the useful spectral region --------

fmin = fc - max(4 * fm, 5000)
fmax = fc + max(4 * fm, 5000)

mask = (freq >= fmin) & (freq <= fmax)

# Convert Vpk spectrum to dBV peak
floor = 1e-12

db_in = 20 * np.log10(
    np.maximum(amp_in[mask], floor)
)

db_out = 20 * np.log10(
    np.maximum(amp_out[mask], floor)
)

plt.figure(figsize=(12, 6))

plt.plot(
    freq[mask] / 1000,
    db_in,
    label="CH1 - Vin"
)

plt.plot(
    freq[mask] / 1000,
    db_out,
    label="CH2 - Vout"
)

plt.axvline(fl / 1000, linestyle="--", alpha=0.5)
plt.axvline(fc / 1000, linestyle="--", alpha=0.5)
plt.axvline(fu / 1000, linestyle="--", alpha=0.5)

plt.xlabel("Frequency (kHz)")
plt.ylabel("Amplitude (dBV peak)")
plt.title(
    f"FFT around AM carrier: "
    f"fc={fc/1000:.3f} kHz, "
    f"fm={fm/1000:.3f} kHz"
)

plt.legend()
plt.grid(True, alpha=0.25)
plt.tight_layout()
plt.show()


scope.close()
gen.close()
rm.close()
