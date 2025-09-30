import pandas as pd
import matplotlib.pyplot as plt

def plot_sequence(file,col,title):
    data = pd.read_csv(file)
    plt.figure(figsize=(6,3))
    plt.stem(data['n'], data[data.columns[1]])
    plt.title(title)
    plt.xlabel('n')
    plt.ylabel(data.columns[1])
    plt.grid(True)

# FIR filters
plot_sequence('fir_lowpass.csv','fir_lowpass','FIR Low-pass')
plot_sequence('fir_highpass.csv','fir_highpass','FIR High-pass')
plot_sequence('fir_bandpass.csv','fir_bandpass','FIR Band-pass')

# Input/output signals
plot_sequence('input_signal.csv','value','Input Signal')
plot_sequence('output_lp.csv','value','Output LP FIR')
plot_sequence('output_hp.csv','value','Output HP FIR')
plot_sequence('output_bp.csv','value','Output BP FIR')
plot_sequence('iir_first_order.csv','value','First-order IIR LP')

plt.show()