import pandas as pd
import matplotlib.pyplot as plt

def plot_sequence(file,title):
    data = pd.read_csv(file)
    plt.figure(figsize=(6,3))
    plt.plot(data['n'], data.iloc[:,1])
    plt.title(title)
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.grid(True)

plot_sequence('uniform_sequence.csv','Uniform Random Sequence')
plot_sequence('gaussian_sequence.csv','Gaussian Random Sequence')
plot_sequence('filtered_gaussian.csv','Filtered Gaussian Sequence')
plot_sequence('autocorrelation.csv','Autocorrelation')
plot_sequence('crosscorrelation.csv','Cross-correlation with Deterministic Sequence')
plot_sequence('exponential_average.csv','Exponential Average of Gaussian Sequence')

plt.show()