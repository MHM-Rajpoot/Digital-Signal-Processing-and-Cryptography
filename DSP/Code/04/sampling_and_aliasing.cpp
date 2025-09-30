/* 
 * sampling_and_aliasing.cpp
 *
 * Demonstration of discrete sequences and spectra:
 *  - Discrete Fourier transform of sequences
 *  - Periodic sampling of continuous signals
 *  - Aliasing due to under-sampling
 *  - Sampling and reconstruction of low-pass and band-pass signals
 *  - Spectral inversion
 *
 * Exports CSV files for plotting in MATLAB, Python, or GNUplot.
 *
 * Written in the style of Computer Laboratory examples
 * (inspired by Dr Markus Kuhn's teaching material).
 *
 * Compile with: g++ -std=c++11 -o sampling_and_aliasing sampling_and_aliasing.cpp
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

/* Utility: export vector<complex<double>> to CSV */
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

/* Utility: export vector<double> to CSV */
void export_csv(const string& filename, const vector<double>& data) {
    ofstream file(filename);
    file << "n,value\n";
    for (size_t n = 0; n < data.size(); ++n) {
        file << n << "," << data[n] << "\n";
    }
    file.close();
    cout << "Exported " << filename << endl;
}

/* Discrete Fourier Transform */
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

/* Inverse DFT */
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

/* Continuous-time sinusoid sampled at fs */
vector<double> sample_sinusoid(double freq, double fs, int N) {
    vector<double> x(N, 0);
    double Ts = 1.0 / fs;
    for (int n = 0; n < N; ++n) {
        x[n] = cos(2.0 * M_PI * freq * n * Ts);
    }
    return x;
}

/* Low-pass reconstruction (ideal sinc interpolation approximation) */
vector<double> reconstruct_lowpass(const vector<double>& x, double fs, int factor) {
    vector<double> y;
    for (size_t n = 0; n < x.size(); ++n) {
        for (int k = 0; k < factor; ++k) {
            y.push_back(x[n]); // simple ZOH upsampling
        }
    }
    return y;
}

/* Spectral inversion: multiply time-domain signal by (-1)^n */
vector<double> spectral_inversion(const vector<double>& x) {
    vector<double> y(x.size(), 0.0);
    for (size_t n = 0; n < x.size(); ++n) {
        y[n] = ((n % 2 == 0) ? 1.0 : -1.0) * x[n];
    }
    return y;
}

int main() {
    cout << "=== Discrete Sequences and Spectra Demonstration ===\n" << endl;

    int N = 32;       // number of samples
    double fs = 16.0; // sampling frequency (Hz)

    /* 1. Discrete sequence and spectrum */
    vector<double> seq = sample_sinusoid(2.0, fs, N); // 2 Hz cosine
    vector<complex<double>> seq_c(seq.begin(), seq.end());
    auto X = DFT(seq_c);
    export_csv("sequence.csv", seq);
    export_csv("sequence_spectrum.csv", X);

    /* 2. Aliasing demonstration */
    auto seq_alias = sample_sinusoid(10.0, fs, N); // 10 Hz cosine (aliased)
    vector<complex<double>> seq_alias_c(seq_alias.begin(), seq_alias.end());
    auto X_alias = DFT(seq_alias_c);
    export_csv("aliased_sequence.csv", seq_alias);
    export_csv("aliased_spectrum.csv", X_alias);

    /* 3. Low-pass reconstruction (upsampling) */
    auto recon_low = reconstruct_lowpass(seq, fs, 4);
    export_csv("lowpass_reconstruction.csv", recon_low);

    /* 4. Band-pass sampling */
    auto seq_band = sample_sinusoid(6.0, fs, N); // 6 Hz (within Nyquist)
    vector<complex<double>> seq_band_c(seq_band.begin(), seq_band.end());
    auto X_band = DFT(seq_band_c);
    export_csv("bandpass_sequence.csv", seq_band);
    export_csv("bandpass_spectrum.csv", X_band);

    /* 5. Spectral inversion */
    auto inverted = spectral_inversion(seq);
    vector<complex<double>> inverted_c(inverted.begin(), inverted.end());
    auto X_inv = DFT(inverted_c);
    export_csv("spectral_inversion.csv", inverted);
    export_csv("spectral_inversion_spectrum.csv", X_inv);

    cout << "\nAll CSV files exported. Plot them to observe spectra, aliasing, and reconstruction." << endl;
    return 0;
}

/* ================================================================
   Python plotting (save as plot_sampling_aliasing.py)
   ================================================================

import pandas as pd
import matplotlib.pyplot as plt

def plot_sequence(csvfile, title):
    df = pd.read_csv(csvfile)
    plt.figure()
    plt.stem(df['n'], df['value'], use_line_collection=True)
    plt.title(title)
    plt.xlabel("n")
    plt.ylabel("x[n]")
    plt.grid(True)

def plot_spectrum(csvfile, title):
    df = pd.read_csv(csvfile)
    plt.figure()
    plt.stem(df['n'], df['mag'], use_line_collection=True)
    plt.title(title)
    plt.xlabel("Frequency bin (k)")
    plt.ylabel("|X[k]|")
    plt.grid(True)

# Original sequence and spectrum
plot_sequence("sequence.csv", "Original 2 Hz Sequence")
plot_spectrum("sequence_spectrum.csv", "Spectrum of 2 Hz Sequence")

# Aliasing demonstration
plot_sequence("aliased_sequence.csv", "Aliased 10 Hz Sequence (fs=16 Hz)")
plot_spectrum("aliased_spectrum.csv", "Spectrum of Aliased 10 Hz Sequence")

# Low-pass reconstruction
plot_sequence("lowpass_reconstruction.csv", "Low-pass Reconstruction (Upsampled)")

# Band-pass sampling
plot_sequence("bandpass_sequence.csv", "Band-pass Sampled 6 Hz Sequence")
plot_spectrum("bandpass_spectrum.csv", "Spectrum of 6 Hz Sequence")

# Spectral inversion
plot_sequence("spectral_inversion.csv", "Spectral Inversion Sequence")
plot_spectrum("spectral_inversion_spectrum.csv", "Spectrum after Spectral Inversion")

plt.show()

================================================================ */
