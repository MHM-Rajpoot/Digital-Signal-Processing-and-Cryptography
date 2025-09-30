import pandas as pd
import matplotlib.pyplot as plt

def plot_sequence(csvfile, title):
    df = pd.read_csv(csvfile)
    plt.figure()
    plt.stem(df['n'], df['value'])
    plt.title(title)
    plt.xlabel("n")
    plt.ylabel("x[n]")
    plt.grid(True)

def plot_spectrum(csvfile, title):
    df = pd.read_csv(csvfile)
    plt.figure()
    plt.stem(df['n'], df['mag'])
    plt.title(title)
    plt.xlabel("Frequency bin (k)")
    plt.ylabel("|X[k]|")
    plt.grid(True)

# Original sequence and spectrum
plot_sequence("sequence.csv", "Original 2 Hz Sequence")
plot_spectrum("sequence_spectrum.csv", "Spectrum of 2 Hz Sequence")

# Aliasing demonstration
plot_sequence("aliased_sequence.csv", "Aliased 10 Hz Sequence (fs=16 Hz)")
plot_spectrum("aliased_spectrum.csv", "Spectrum of Aliased 10 Hz Sequence")

# Low-pass reconstruction
plot_sequence("lowpass_reconstruction.csv", "Low-pass Reconstruction (Upsampled)")

# Band-pass sampling
plot_sequence("bandpass_sequence.csv", "Band-pass Sampled 6 Hz Sequence")
plot_spectrum("bandpass_spectrum.csv", "Spectrum of 6 Hz Sequence")

# Spectral inversion
plot_sequence("spectral_inversion.csv", "Spectral Inversion Sequence")
plot_spectrum("spectral_inversion_spectrum.csv", "Spectrum after Spectral Inversion")

plt.show()
