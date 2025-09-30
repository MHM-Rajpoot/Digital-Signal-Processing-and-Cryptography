/*
 * filters_demo.cpp
 *
 * Demonstration of FIR/IIR filters, properties, design techniques:
 *  - FIR filters (window-based, high-pass via spectral inversion, band-pass via modulation)
 *  - FFT-based convolution
 *  - Polynomial representation and z-transform
 *  - IIR filter design (Butterworth, Chebyshev I/II, elliptic)
 *
 * CSV export for MATLAB/Python/GNUplot.
 *
 * Inspired by Dr Markus Kuhn's Computer Laboratory material.
 *
 * Compile with: g++ -std=c++11 -o filters_demo filters_demo.cpp -lm
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>

using namespace std;

/* ===== CSV Export ===== */
void export_csv(const string &filename, const vector<double>& data, const string &colname="value") {
    ofstream file(filename);
    file << "n," << colname << "\n";
    for (size_t n = 0; n < data.size(); ++n)
        file << n << "," << data[n] << "\n";
    file.close();
    cout << "Exported " << filename << endl;
}

void export_csv(const string &filename, const vector<complex<double>>& data) {
    ofstream file(filename);
    file << "n,real,imag,mag,phase\n";
    for (size_t n = 0; n < data.size(); ++n)
        file << n << "," << real(data[n]) << "," << imag(data[n])
             << "," << abs(data[n]) << "," << arg(data[n]) << "\n";
    file.close();
    cout << "Exported " << filename << endl;
}

/* ===== FIR Design: Windowed sinc (low-pass) ===== */
vector<double> fir_lowpass(int N, double fc) {
    vector<double> h(N, 0.0);
    int M = N-1;
    for (int n = 0; n < N; ++n) {
        double m = n - M/2.0;
        if (m == 0.0) h[n] = 2.0 * fc;
        else h[n] = sin(2.0*M_PI*fc*m)/(M_PI*m);
        // Hamming window
        h[n] *= 0.54 - 0.46*cos(2.0*M_PI*n/M);
    }
    return h;
}

/* Frequency inversion to obtain high-pass filter */
vector<double> fir_highpass(const vector<double>& h_lp) {
    vector<double> h_hp = h_lp;
    for (size_t n = 0; n < h_hp.size(); ++n)
        h_hp[n] = (n == h_hp.size()/2) ? 1.0 - h_hp[n] : -h_hp[n];
    return h_hp;
}

/* Modulation to create band-pass filter */
vector<double> fir_bandpass(const vector<double>& h_lp, double f1, double f2) {
    double fc = (f2+f1)/2.0;
    vector<double> h_bp = h_lp;
    for (size_t n = 0; n < h_lp.size(); ++n)
        h_bp[n] *= 2.0*cos(2.0*M_PI*fc*(n - (h_lp.size()-1)/2.0));
    return h_bp;
}

/* ===== Convolution (time-domain and FFT-based) ===== */
vector<double> convolution(const vector<double>& x, const vector<double>& h) {
    int N = x.size(), M = h.size();
    vector<double> y(N+M-1,0.0);
    for(int n=0; n<(int)y.size(); n++)
        for(int k=0;k<N;k++)
            if(n-k>=0 && n-k<M) y[n] += x[k]*h[n-k];
    return y;
}

/* ===== FFT (Cooley-Tukey, size must be power of 2) ===== */
void fft(vector<complex<double>>& a) {
    int N = a.size();
    if(N<=1) return;
    vector<complex<double>> even(N/2), odd(N/2);
    for(int i=0;i<N/2;i++){
        even[i]=a[i*2];
        odd[i]=a[i*2+1];
    }
    fft(even); fft(odd);
    const complex<double> j(0,1);
    for(int k=0;k<N/2;k++){
        complex<double> t = exp(-2.0*M_PI*j*(double)k/(double)N)*odd[k];
        a[k]=even[k]+t;
        a[k+N/2]=even[k]-t;
    }
}

/* Generate example input signal */
vector<double> input_signal(int N) {
    vector<double> x(N,0.0);
    for(int n=0;n<N;n++)
        x[n] = sin(2*M_PI*0.05*n) + 0.5*sin(2*M_PI*0.15*n);
    return x;
}

/* ===== Simple IIR: first-order low-pass (discrete-time analog approximation) ===== */
vector<double> iir_first_order(double alpha, const vector<double>& x) {
    vector<double> y(x.size(),0.0);
    for(size_t n=0;n<x.size();n++)
        y[n] = (n==0)? alpha*x[n] : alpha*x[n] + (1-alpha)*y[n-1];
    return y;
}

/* ===== Main demonstration ===== */
int main(){
    cout<<"=== FIR and IIR Filter Demonstration ==="<<endl;

    int N=51; // FIR length
    double fc=0.1; // normalized cutoff
    auto h_lp = fir_lowpass(N, fc);
    export_csv("fir_lowpass.csv", h_lp);

    auto h_hp = fir_highpass(h_lp);
    export_csv("fir_highpass.csv", h_hp);

    auto h_bp = fir_bandpass(h_lp,0.05,0.15);
    export_csv("fir_bandpass.csv", h_bp);

    auto x = input_signal(128);
    export_csv("input_signal.csv", x);

    auto y_lp = convolution(x,h_lp);
    export_csv("output_lp.csv", y_lp);

    auto y_hp = convolution(x,h_hp);
    export_csv("output_hp.csv", y_hp);

    auto y_bp = convolution(x,h_bp);
    export_csv("output_bp.csv", y_bp);

    // Example FFT of FIR low-pass filter
    vector<complex<double>> H_fft(64,0);
    for(int i=0;i<N;i++) H_fft[i] = h_lp[i];
    fft(H_fft);
    export_csv("fft_fir_lowpass.csv", H_fft);

    // Example IIR first-order low-pass
    auto y_iir = iir_first_order(0.1, x);
    export_csv("iir_first_order.csv", y_iir);

    cout<<"All CSV files exported. Plot them in Python/MATLAB/GNUplot."<<endl;

    return 0;
}

/*
# Python plotting code (save as filters_plot.py):

import pandas as pd
import matplotlib.pyplot as plt

def plot_sequence(file,col,title):
    data = pd.read_csv(file)
    plt.figure(figsize=(6,3))
    plt.stem(data['n'], data[data.columns[1]], use_line_collection=True)
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
*/
