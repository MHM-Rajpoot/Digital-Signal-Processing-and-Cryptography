import pandas as pd
import matplotlib.pyplot as plt

def plot_complex_csv(file,title):
    data = pd.read_csv(file)
    plt.figure(figsize=(6,3))
    plt.plot(data['n'], data['real'], label='Real')
    plt.plot(data['n'], data['imag'], label='Imag')
    plt.title(title)
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.legend()
    plt.grid(True)

plot_complex_csv('am_signal.csv','AM Signal')
plot_complex_csv('fm_signal.csv','FM Signal')
plot_complex_csv('msk_signal.csv','MSK Signal')
plot_complex_csv('qam_signal.csv','QAM Signal')
plot_complex_csv('ofdm_signal.csv','OFDM Signal')

# Matched filter output
mf = pd.read_csv('matched_filter_output.csv')
plt.figure(figsize=(6,3))
plt.stem(mf['n'], mf['value'])
plt.title('Matched Filter Output')
plt.xlabel('n')
plt.ylabel('Correlation Magnitude')
plt.grid(True)

plt.show()