import pandas as pd
import matplotlib.pyplot as plt

def plot_sequence(filename, colname, title):
    data = pd.read_csv(filename)
    plt.figure(figsize=(6,4))
    plt.stem(data['n'], data[colname], basefmt='k-')
    plt.title(title)
    plt.xlabel('n')
    plt.ylabel(colname)
    plt.grid(True)

# Plot unit impulse δ[n]
plot_sequence('delta.csv', 'δ[n]', 'Unit Impulse δ[n]')

# Plot unit step u[n-2]
plot_sequence('step.csv', 'u[n-2]', 'Unit Step u[n-2]')

# Plot x[n], h[n], and y[n] from convolution
plot_sequence('x.csv', 'x[n]', 'Signal x[n]')
plot_sequence('h.csv', 'h[n]', 'Impulse Response h[n]')
plot_sequence('y_conv.csv', 'y[n]', 'Convolution y[n] = x[n] * h[n]')

# Plot moving average system input/output
plot_sequence('input.csv', 'input[n]', 'Input Signal to Moving Average')
plot_sequence('moving_avg.csv', 'y[n] (moving average)', 'Moving Average Output (LTI System)')

plt.tight_layout()
plt.show()