import pandas as pd
import matplotlib.pyplot as plt

orig = pd.read_csv('original_signal.csv')['value']
quant = pd.read_csv('quantized_signal.csv')['value']
masked = pd.read_csv('masked_signal.csv')['value']

plt.figure(figsize=(12,6))
plt.plot(orig, label='Original Signal')
plt.plot(quant, label='Lossy Quantized', linestyle='--')
plt.plot(masked, label='Masked/Reconstructed', linestyle=':')
plt.title("Lossless vs Lossy Compression and Psychoacoustic Masking")
plt.xlabel("Sample Index")
plt.ylabel("Amplitude")
plt.grid(True)
plt.legend()
plt.show()