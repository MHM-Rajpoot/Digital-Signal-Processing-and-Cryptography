import pandas as pd
import matplotlib.pyplot as plt

def plot_spectrum(csvfile, title):
    df = pd.read_csv(csvfile)
    plt.figure()
    plt.stem(df['n'], df['mag'])
    plt.title(title)
    plt.xlabel("Frequency bin (n)")
    plt.ylabel("|X[n]|")
    plt.grid(True)

# Cosine spectrum
plot_spectrum("cosine_spectrum.csv", "Cosine Spectrum (DFT)")

# Linearity check
plot_spectrum("linearity_sum_spectrum.csv", "Linearity: Spectrum of (x1+x2)")
plot_spectrum("linearity_expected_spectrum.csv", "Linearity: Spectrum of x1 + Spectrum of x2")

# FFT vs DFT
plot_spectrum("fft_spectrum.csv", "FFT Spectrum")

# Real FFT
plot_spectrum("real_fft_spectrum.csv", "Real-valued FFT Spectrum")

plt.show()
