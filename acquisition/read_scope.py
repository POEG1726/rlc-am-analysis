
# Local instrument configuration is required before legacy execution.
from acquisition_config import require_context
_bench = require_context()

import pyvisa
import numpy as np
import struct
import matplotlib.pyplot as plt
import pandas as pd

SCOPE_ADDR = _bench.scope_address

rm = pyvisa.ResourceManager()
scope = rm.open_resource(SCOPE_ADDR)

scope.timeout = 15000
scope.chunk_size = 1024 * 1024


def query_binary_block(inst, command):
    """
    Send a SCPI query returning an IEEE488.2 binary block:
    #<n><length><binary data>
    """
    inst.write(command)
    raw = inst.read_raw()

    # SDS800X HD can occasionally prepend text before '#'
    pos = raw.find(b"#")
    if pos < 0:
        raise RuntimeError(
            f"No binary block header found for {command!r}. "
            f"First bytes: {raw[:100]!r}"
        )

    n_digits = int(chr(raw[pos + 1]))

    length_start = pos + 2
    length_end = length_start + n_digits

    data_length = int(raw[length_start:length_end].decode())
    data_start = length_end
    data_end = data_start + data_length

    return raw[data_start:data_end]


def read_channel(ch, max_points=200_000):
    # Select waveform source
    scope.write(f":WAV:SOUR C{ch}")

    # SDS824X HD = 12-bit ADC, therefore use WORD
    scope.write(":WAV:WIDT WORD")
    scope.write(":WAV:BYT LSB")

    # Start from first point, no decimation
    scope.write(":WAV:STAR 0")
    scope.write(":WAV:INT 1")

    # Don't ask for more points than currently acquired
    acq_points = int(float(scope.query(":ACQ:POIN?")))
    points = min(acq_points, max_points)

    scope.write(f":WAV:POIN {points}")

    # Get waveform descriptor / preamble
    pre = query_binary_block(scope, ":WAV:PRE?")

    if len(pre) < 188:
        raise RuntimeError(
            f"Preamble unexpectedly short: {len(pre)} bytes"
        )

    # SDS waveform descriptor fields
    vdiv = struct.unpack_from("<f", pre, 156)[0]
    voffset = struct.unpack_from("<f", pre, 160)[0]
    code_per_div = struct.unpack_from("<f", pre, 164)[0]
    dt = struct.unpack_from("<f", pre, 176)[0]

    # Get waveform samples
    data = query_binary_block(scope, ":WAV:DATA?")

    # 16-bit signed, little-endian
    codes = np.frombuffer(data, dtype="<i2")

    # Convert ADC codes to volts
    volts = (
        codes.astype(np.float64)
        * (vdiv / code_per_div)
        - voffset
    )

    # Time axis
    t = np.arange(len(volts)) * dt

    info = {
        "channel": ch,
        "points": len(volts),
        "dt": dt,
        "sample_rate": 1.0 / dt,
        "vdiv": vdiv,
        "voffset": voffset,
        "code_per_div": code_per_div,
    }

    return t, volts, info


print("Scope:", scope.query("*IDN?").strip())

print("Acquisition points:", scope.query(":ACQ:POIN?").strip())
print("Sample rate:", scope.query(":ACQ:SRAT?").strip())
print("Trigger state:", scope.query(":TRIG:STAT?").strip())

# Freeze acquisition so CH1 and CH2 are read from the same captured record
scope.write(":TRIG:STOP")

print("\nReading CH1...")
t1, v1, info1 = read_channel(1)

print("Reading CH2...")
t2, v2, info2 = read_channel(2)

print("\nCH1:")
for k, v in info1.items():
    print(f"  {k}: {v}")

print("\nCH2:")
for k, v in info2.items():
    print(f"  {k}: {v}")

# Resume scope
scope.write(":TRIG:RUN")

# Use common number of points
n = min(len(v1), len(v2))

df = pd.DataFrame({
    "time_s": t1[:n],
    "vin_ch1_V": v1[:n],
    "vout_ch2_V": v2[:n],
})

df.to_csv("scope_test.csv", index=False)

print(f"\nSaved {n} points to scope_test.csv")

# Quick visual check
plt.figure(figsize=(12, 6))
plt.plot(t1[:n] * 1000, v1[:n], label="CH1 - Vin", linewidth=0.8)
plt.plot(t2[:n] * 1000, v2[:n], label="CH2 - Vout", linewidth=0.8)

plt.xlabel("Time (ms)")
plt.ylabel("Voltage (V)")
plt.legend()
plt.tight_layout()
plt.show()

scope.close()
rm.close()
