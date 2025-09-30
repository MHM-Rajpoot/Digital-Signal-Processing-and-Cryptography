/* 
 * fourier_transform_demo.cpp
 *
 * Demonstration of Fourier transform, phasors as orthogonal basis functions,
 * forms of the Fourier transform, convolution theorem, Dirac’s delta,
 * and impulse combs in time/frequency domain.
 *
 * Includes CSV export for plotting in MATLAB, Python, or GNUplot.
 *
 * Written in the style of Computer Laboratory examples
 * (inspired by Dr Markus Kuhn's teaching material).
 *
 * Compile with: g++ -std=c++11 -o fourier_transform_demo fourier_transform_demo.cpp
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

/* Utility: export vector to CSV */
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

/* Example 1: Discrete Fourier Transform (DFT) */
vector<complex<double>> DFT(const vector<complex<double>>& x) {
    int N = x.size();
    vector<complex<double>> X(N, 0);
    const complex<double> j(0, 1);

    for (int k = 0; k < N; ++k) {
        for (int n = 0; n < N; ++n) {
            double angle = -2 * M_PI * k * n / N;
            X[k] += x[n] * exp(j * angle);
        }
    }
    return X;
}

/* Example 2: Inverse DFT */
vector<complex<double>> IDFT(const vector<complex<double>>& X) {
    int N = X.size();
    vector<complex<double>> x(N, 0);
    const complex<double> j(0, 1);

    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < N; ++k) {
            double angle = 2 * M_PI * k * n / N;
            x[n] += X[k] * exp(j * angle);
        }
        x[n] /= (double)N;
    }
    return x;
}

/* Example 3: Convolution theorem demonstration */
vector<complex<double>> convolution(const vector<complex<double>>& x, const vector<complex<double>>& h) {
    int N = x.size();
    int M = h.size();
    vector<complex<double>> y(N + M - 1, 0);

    for (int n = 0; n < (int)y.size(); ++n) {
        for (int k = 0; k < N; ++k) {
            if (n - k >= 0 && n - k < M)
                y[n] += x[k] * h[n - k];
        }
    }
    return y;
}

/* Example 4: Dirac delta approximation */
vector<complex<double>> dirac_delta(int N, int pos) {
    vector<complex<double>> delta(N, 0);
    if (pos >= 0 && pos < N)
        delta[pos] = 1.0;
    return delta;
}

/* Example 5: Impulse train (comb) */
vector<complex<double>> impulse_train(int N, int period) {
    vector<complex<double>> comb(N, 0);
    for (int n = 0; n < N; n += period) {
        comb[n] = 1.0;
    }
    return comb;
}

int main() {
    cout << "=== Fourier Transform Demonstration ===\n" << endl;

    /* Fourier transform of a simple sequence */
    vector<complex<double>> x = {1, 1, 1, 0, 0, 0, 0, 0}; // Rectangular pulse
    auto X = DFT(x);
    export_csv("rectangular_pulse.csv", x);
    export_csv("rectangular_pulse_spectrum.csv", X);

    /* Reconstruction with IDFT */
    auto x_recon = IDFT(X);
    export_csv("rectangular_pulse_reconstructed.csv", x_recon);

    /* Phasors as orthogonal basis functions */
    cout << "\n=== Phasors as Orthogonal Basis Functions ===" << endl;
    const complex<double> j(0, 1);
    int N = 8;
    complex<double> inner_prod = 0;
    for (int n = 0; n < N; ++n) {
        double n_over_N = static_cast<double>(n) / static_cast<double>(N);
        inner_prod += exp(j * 2.0 * M_PI * 1.0 * n_over_N) *
                    conj(exp(j * 2.0 * M_PI * 2.0 * n_over_N));
    }
    cout << "Inner product of e^(j2πn/N) and e^(j4πn/N) over N=8 = "
        << inner_prod << " (≈0 → orthogonal)" << endl;

    /* Convolution theorem */
    cout << "\n=== Convolution Theorem ===" << endl;
    vector<complex<double>> h = {1, 0.5, 0.25};
    auto y_time = convolution(x, h);
    auto X1 = DFT(x);
    auto H = DFT(h);
    vector<complex<double>> Y_freq(X1.size());
    for (size_t k = 0; k < X1.size(); ++k)
        Y_freq[k] = X1[k] * H[k];
    auto y_freq = IDFT(Y_freq);

    export_csv("convolution_time.csv", y_time);
    export_csv("convolution_freq.csv", y_freq);

    /* Dirac delta */
    auto delta = dirac_delta(8, 0);
    export_csv("dirac_delta.csv", delta);
    auto D = DFT(delta);
    export_csv("dirac_delta_spectrum.csv", D);

    /* Impulse train */
    auto comb = impulse_train(16, 4);
    export_csv("impulse_train.csv", comb);
    auto C = DFT(comb);
    export_csv("impulse_train_spectrum.csv", C);
    
    cout << "\nAll CSV files exported. Use MATLAB, Python, or GNUplot for plotting.\n";
    return 0;
}

/* =====================================================================
   Python plotting code (save separately as plot_fourier_demo.py)
   =====================================================================

import pandas as pd
import matplotlib.pyplot as plt

def plot_csv(filename, title, time_domain=True):
    data = pd.read_csv(filename)
    plt.figure(figsize=(8,4))
    if time_domain:
        plt.stem(data['n'], data['real'], use_line_collection=True)
        plt.title(title + " (Time Domain)")
        plt.xlabel("n")
        plt.ylabel("Value")
    else:
        plt.stem(data['n'], data['mag'], use_line_collection=True)
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

===================================================================== */
