import pandas as pd
import matplotlib.pyplot as plt

def plot_csv(filename, title, time_domain=True):
    data = pd.read_csv(filename)
    plt.figure(figsize=(8,4))
    if time_domain:
        plt.stem(data['n'], data['real'])
        plt.title(title + " (Time Domain)")
        plt.xlabel("n")
        plt.ylabel("Value")
    else:
        plt.stem(data['n'], data['mag'])
        plt.title(title + " (Spectrum Magnitude)")
        plt.xlabel("Frequency bin k")
        plt.ylabel("|X[k]|")
    plt.grid(True)

# Plot examples
plot_csv("rectangular_pulse.csv", "Rectangular Pulse", time_domain=True)
plot_csv("rectangular_pulse_spectrum.csv", "Rectangular Pulse Spectrum", time_domain=False)

plot_csv("convolution_time.csv", "Convolution in Time", time_domain=True)
plot_csv("convolution_freq.csv", "Convolution via Frequency Domain", time_domain=True)

plot_csv("dirac_delta.csv", "Dirac Delta", time_domain=True)
plot_csv("dirac_delta_spectrum.csv", "Dirac Delta Spectrum", time_domain=False)

plot_csv("impulse_train.csv", "Impulse Train", time_domain=True)
plot_csv("impulse_train_spectrum.csv", "Impulse Train Spectrum", time_domain=False)

plt.show()
