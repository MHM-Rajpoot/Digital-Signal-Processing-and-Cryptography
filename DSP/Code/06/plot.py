import pandas as pd
import matplotlib.pyplot as plt

def plot_spectrum(csvfile, title):
    df = pd.read_csv(csvfile)
    plt.figure()
    plt.stem(df['k'], df['mag'])
    plt.title(title)
    plt.xlabel("Frequency bin (k)")
    plt.ylabel("|X[k]|")
    plt.grid(True)

# Leakage
plot_spectrum("leakage_spectrum.csv", "Leakage (f=5.3 Hz)")

# Scalloping loss
plot_spectrum("scalloping_f5_spectrum.csv", "Scalloping (f=5 Hz)")
plot_spectrum("scalloping_f5p5_spectrum.csv", "Scalloping (f=5.5 Hz)")

# Windowing
plot_spectrum("window_rectangular_spectrum.csv", "Rectangular Window")
plot_spectrum("window_hann_spectrum.csv", "Hann Window")

# Zero padding
plot_spectrum("zero_padding_spectrum.csv", "Zero Padding (N=256)")

plt.show()
