import pandas as pd
import matplotlib.pyplot as plt

def plot_time_signal(filename, title):
    data = pd.read_csv(filename)
    plt.figure(figsize=(6,4))
    plt.stem(data['n'], data['real'], linefmt='b-', markerfmt='bo', basefmt='k-')
    plt.title(title + " (Real Part)")
    plt.xlabel("n")
    plt.ylabel("Value")
    plt.grid(True)

def plot_phasors(filename, title):
    data = pd.read_csv(filename)
    plt.figure(figsize=(6,6))
    plt.quiver([0]*len(data), [0]*len(data),
               data['real'], data['imag'],
               angles='xy', scale_units='xy', scale=1, color=['r','g','b'])
    for i, txt in enumerate(data.index):
        plt.text(data['real'][i]*1.05, data['imag'][i]*1.05, f"{i}")
    plt.title(title + " (Phasors)")
    plt.xlabel("Re")
    plt.ylabel("Im")
    plt.grid(True)
    plt.axis('equal')

# Plotting
plot_phasors("phasors.csv", "Phasor Demo")
plot_time_signal("eigen_x.csv", "Eigenfunction Input x[n]")
plot_time_signal("eigen_y.csv", "Eigenfunction Output y[n]")
plot_phasors("applications.csv", "Applications (Electronics, Optics, Acoustics)")

plt.show()
