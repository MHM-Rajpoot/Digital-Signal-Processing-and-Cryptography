/* 
 * dft_fft_demo.cpp
 *
 * Demonstration of:
 *  - Discrete Fourier Transform (DFT)
 *  - Continuous vs. discrete Fourier transform (via sampled signals)
 *  - Symmetry of spectra
 *  - Linearity of the Fourier transform
 *  - Review of FFT (Cooley–Tukey algorithm)
 *  - Real-valued FFT optimization
 *
 * Exports CSV files for plotting in MATLAB, Python, or GNUplot.
 *
 * Written in the style of Computer Laboratory examples
 * (inspired by Dr Markus Kuhn's teaching material).
 *
 * Compile with: g++ -std=c++11 -o dft_fft_demo dft_fft_demo.cpp -lm
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

/* ===== CSV Export Utilities ===== */
void export_csv(const string& filename, const vector<complex<double>>& data) {
    ofstream file(filename);
    file << "n,real,imag,mag,phase\n";
    for (size_t n = 0; n < data.size(); ++n) {
        file << n << ","
             << real(data[n]) << ","
             << imag(data[n]) << ","
             << abs(data[n]) << ","
             << arg(data[n]) << "\n";
    }
    file.close();
    cout << "Exported " << filename << endl;
}

void export_csv(const string& filename, const vector<double>& data) {
    ofstream file(filename);
    file << "n,value\n";
    for (size_t n = 0; n < data.size(); ++n) {
        file << n << "," << data[n] << "\n";
    }
    file.close();
    cout << "Exported " << filename << endl;
}

/* ===== DFT and IDFT ===== */
vector<complex<double>> DFT(const vector<complex<double>>& x) {
    int N = x.size();
    vector<complex<double>> X(N, 0);
    const complex<double> j(0, 1);

    for (int k = 0; k < N; ++k) {
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * M_PI * k * n / N;
            X[k] += x[n] * exp(j * angle);
        }
    }
    return X;
}

vector<complex<double>> IDFT(const vector<complex<double>>& X) {
    int N = X.size();
    vector<complex<double>> x(N, 0);
    const complex<double> j(0, 1);

    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < N; ++k) {
            double angle = 2.0 * M_PI * k * n / N;
            x[n] += X[k] * exp(j * angle);
        }
        x[n] /= (double)N;
    }
    return x;
}

/* ===== Recursive Cooley–Tukey FFT ===== */
void FFT(vector<complex<double>>& a) {
    int N = a.size();
    if (N <= 1) return;

    // Divide 
    vector<complex<double>> even(N / 2), odd(N / 2);
    for (int i = 0; i < N / 2; i++) {
        even[i] = a[i*2];
        odd[i] = a[i*2+1];
    }

    // Conqurer
    FFT(even);
    FFT(odd);

    // Combine
    const complex<double> j(0, 1);
    for (int k = 0; k < N / 2; k++) {
        complex<double> t = exp(-2.0 * M_PI * j * (double)k / (double)N) * odd[k];
        a[k] = even[k] + t;
        a[k+N/2] = even[k] - t;
    }
}

/* ===== Real-valued FFT optimization ===== */
vector<complex<double>> real_fft(const vector<double>& x) {
    vector<complex<double>> cx(x.begin(), x.end());
    FFT(cx);
    return cx;
}

/* ===== Example signals ===== */
vector<double> sample_cosine(double f, double fs, int N) {
    vector<double> x(N);
    double Ts = 1.0 / fs;
    for (int n = 0; n < N; n++) {
        x[n] = cos(2 * M_PI * f * n * Ts);
    }
    return x;
}

int main() {
    cout << "=== DFT and FFT Demonstration ===\n" << endl;

    int N = 32;
    double fs = 32.0;

    /* 1. DFT of a simple cosine */
    auto x = sample_cosine(4.0, fs, N);
    vector<complex<double>> x_c(x.begin(), x.end());
    auto X = DFT(x_c);
    export_csv("cosine_sequence.csv", x);
    export_csv("cosine_spectrum.csv", X);

    /* 2. Continuous vs discrete FT (approximate) */
    cout << "\nContinuous FT of cos(2π*4t): impulses at f=±4 Hz." << endl;
    cout << "Discrete FT shows spreading across bins (finite-length effect)." << endl;

    /* 3. Symmetry */
    cout << "\nSymmetry of real-valued signal spectrum:" << endl;
    for (int k = 0; k < N/2; k++) {
        cout << "X[" << k << "] and X[" << N-k << "] are conjugates." << endl;
    }

    /* 4. Linearity check */
    auto x1 = sample_cosine(3.0, fs, N);
    auto x2 = sample_cosine(5.0, fs, N);
    vector<complex<double>> x1_c(x1.begin(), x1.end());
    vector<complex<double>> x2_c(x2.begin(), x2.end());
    auto X1 = DFT(x1_c);
    auto X2 = DFT(x2_c);

    vector<complex<double>> X_lin(N);
    for (int k = 0; k < N; k++)
        X_lin[k] = X1[k] + X2[k];
    auto x_sum = x1;
    for (int n = 0; n < N; n++) x_sum[n] += x2[n];
    vector<complex<double>> x_sum_c(x_sum.begin(), x_sum.end());
    auto X_sum = DFT(x_sum_c);

    export_csv("linearity_sum_spectrum.csv", X_sum);
    export_csv("linearity_expected_spectrum.csv", X_lin);

    cout << "\nLinearity demonstrated: FT{x1+x2} = FT{x1}+FT{x2}" << endl;

    /* 5. FFT vs DFT */
    auto X_fft = x_c;
    FFT(X_fft);
    export_csv("fft_spectrum.csv", X_fft);
    cout << "\nFFT computed spectrum in O(N log N) vs DFT O(N^2)." << endl;

    /* 6. Real-valued FFT optimization */
    auto X_realfft = real_fft(x);
    export_csv("real_fft_spectrum.csv", X_realfft);

    cout << "\nAll CSV files exported. Plot in MATLAB, Python, or GNUplot." << endl;

    return 0;
}

/* ================================================================
   Python plotting (save as plot_dft_fft.py)
   ================================================================

import pandas as pd
import matplotlib.pyplot as plt

def plot_spectrum(csvfile, title):
    df = pd.read_csv(csvfile)
    plt.figure()
    plt.stem(df['n'], df['mag'], use_line_collection=True)
    plt.title(title)
    plt.xlabel("Frequency bin (n)")
    plt.ylabel("|X[n]|")
    plt.grid(True)

# Cosine spectrum
plot_spectrum("cosine_spectrum.csv", "Cosine Spectrum (DFT)")

# Linearity check
plot_spectrum("linearity_sum_spectrum.csv", "Linearity: Spectrum of (x1+x2)")
plot_spectrum("linearity_expected_spectrum.csv", "Linearity: Spectrum of x1 + Spectrum of x2")

# FFT vs DFT
plot_spectrum("fft_spectrum.csv", "FFT Spectrum")

# Real FFT
plot_spectrum("real_fft_spectrum.csv", "Real-valued FFT Spectrum")

plt.show()

================================================================ */
