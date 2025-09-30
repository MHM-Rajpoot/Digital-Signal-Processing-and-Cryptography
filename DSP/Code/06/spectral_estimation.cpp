/* 
 * spectral_estimation.cpp
 *
 * Demonstration of:
 *   - Spectral leakage
 *   - Scalloping loss
 *   - Windowing (rectangular vs Hann)
 *   - Zero padding
 *
 * Exports CSV files for plotting in MATLAB, Python, or GNUplot.
 *
 * Written in the style of Computer Laboratory examples
 * (inspired by Dr Markus Kuhn's teaching material).
 *
 * Compile with: g++ -std=c++11 -o spectral_estimation spectral_estimation.cpp -lm
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

/* ===== CSV Export ===== */
void export_csv(const string& filename, const vector<complex<double>>& data) {
    ofstream file(filename);
    file << "k,real,imag,mag,phase\n";
    for (size_t k = 0; k < data.size(); ++k) {
        file << k << ","
             << real(data[k]) << ","
             << imag(data[k]) << ","
             << abs(data[k]) << ","
             << arg(data[k]) << "\n";
    }
    file.close();
    cout << "Exported " << filename << endl;
}

void export_csv(const string& filename, const vector<double>& data) {
    ofstream file(filename);
    file << "n,value\n";
    for (size_t n = 0; n < data.size(); ++n)
        file << n << "," << data[n] << "\n";
    file.close();
    cout << "Exported " << filename << endl;
}

/* ===== FFT (recursive Cooley–Tukey) ===== */
void FFT(vector<complex<double>>& a) {
    int N = a.size();
    if (N <= 1) return;

    vector<complex<double>> even(N/2), odd(N/2);
    for (int i = 0; i < N/2; i++) {
        even[i] = a[i*2];
        odd[i]  = a[i*2+1];
    }

    FFT(even);
    FFT(odd);

    const complex<double> j(0,1);
    for (int k = 0; k < N/2; k++) {
        complex<double> t = exp(-2.0*M_PI*j*(double)k/(double)N) * odd[k];
        a[k] = even[k] + t;
        a[k+N/2] = even[k] - t;
    }
}

/* ===== Windows ===== */
vector<double> hann_window(int N) {
    vector<double> w(N);
    for (int n = 0; n < N; n++)
        w[n] = 0.5 * (1 - cos(2*M_PI*n/(N-1)));
    return w;
}

/* ===== Signal generator ===== */
vector<double> sample_cosine(double f, double fs, int N) {
    vector<double> x(N);
    double Ts = 1.0 / fs;
    for (int n = 0; n < N; n++)
        x[n] = cos(2*M_PI*f*n*Ts);
    return x;
}

/* ===== Main ===== */
int main() {
    cout << "=== Spectral Estimation Demonstration ===\n" << endl;

    int N = 64;          // base FFT length
    double fs = 64.0;    // sampling frequency
    double f0 = 5.3;     // sinusoid frequency (not integer bin) → leakage

    /* 1. Leakage demonstration */
    auto x = sample_cosine(f0, fs, N);
    vector<complex<double>> X(x.begin(), x.end());
    FFT(X);
    export_csv("leakage_signal.csv", x);
    export_csv("leakage_spectrum.csv", X);
    cout << "Leakage: sinusoid at " << f0 << " Hz spreads across bins." << endl;

    /* 2. Scalloping loss (compare f=5 and f=5.5 Hz) */
    auto x1 = sample_cosine(5.0, fs, N);
    auto x2 = sample_cosine(5.5, fs, N);
    vector<complex<double>> X1(x1.begin(), x1.end()), X2(x2.begin(), x2.end());
    FFT(X1); FFT(X2);
    export_csv("scalloping_f5_spectrum.csv", X1);
    export_csv("scalloping_f5p5_spectrum.csv", X2);
    cout << "Scalloping: amplitude depends on frequency-bin alignment." << endl;

    /* 3. Windowing (Hann vs Rectangular) */
    auto w = hann_window(N);
    vector<double> x_win(N);
    for (int n = 0; n < N; n++)
        x_win[n] = x[n] * w[n];
    vector<complex<double>> X_rect(x.begin(), x.end());
    vector<complex<double>> X_hann(x_win.begin(), x_win.end());
    FFT(X_rect);
    FFT(X_hann);
    export_csv("window_rectangular_spectrum.csv", X_rect);
    export_csv("window_hann_spectrum.csv", X_hann);
    cout << "Windowing reduces leakage at the cost of resolution." << endl;

    /* 4. Zero padding (increase FFT size to 256) */
    int Npad = 256;
    vector<complex<double>> X_pad(Npad,0);
    for (int n = 0; n < N; n++) X_pad[n] = x[n];
    FFT(X_pad);
    export_csv("zero_padding_spectrum.csv", X_pad);
    cout << "Zero padding: finer spectral grid, no new information." << endl;

    cout << "\nAll CSV files exported. Plot in MATLAB, Python, or GNUplot." << endl;
    return 0;
}

/* ================================================================
   Python plotting (save this as plot_spectral_estimation.py)
   ================================================================

import pandas as pd
import matplotlib.pyplot as plt

def plot_spectrum(csvfile, title):
    df = pd.read_csv(csvfile)
    plt.figure()
    plt.stem(df['k'], df['mag'], use_line_collection=True)
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

================================================================ */
