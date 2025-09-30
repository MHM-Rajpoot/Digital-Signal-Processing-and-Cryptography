/*
 * compression_demo.cpp
 *
 * Demonstration of lossless vs lossy compression and
 * psychoacoustic-inspired lossy compression.
 *
 * Includes CSV export for MATLAB, Python, or GNUplot.
 *
 * Compile with: g++ -std=c++11 -o compression_demo compression_demo.cpp -lm
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;
const double PI = 3.14159265358979323846;

/* ===== CSV export ===== */
void export_csv(const string &filename, const vector<double> &data) {
    ofstream file(filename);
    file << "n,value\n";
    for(size_t i=0;i<data.size();i++) file << i << "," << data[i] << "\n";
    file.close();
    cout << "Exported " << filename << endl;
}

/* ===== Simple DFT/IDFT for psychoacoustic demo ===== */
vector<complex<double>> DFT(const vector<double>& x) {
    int N = x.size();
    vector<complex<double>> X(N);
    for(int k=0;k<N;k++){
        complex<double> sum(0,0);
        for(int n=0;n<N;n++){
            double angle = -2*PI*k*n/N;
            sum += complex<double>(x[n]*cos(angle), x[n]*sin(angle));
        }
        X[k] = sum;
    }
    return X;
}

vector<double> IDFT(const vector<complex<double>>& X) {
    int N = X.size();
    vector<double> x(N,0);
    for(int n=0;n<N;n++){
        complex<double> sum(0,0);
        for(int k=0;k<N;k++){
            double angle = 2*PI*k*n/N;
            sum += X[k]*complex<double>(cos(angle), sin(angle));
        }
        x[n] = sum.real()/N;
    }
    return x;
}

/* ===== Psychoacoustic masking ===== */
vector<complex<double>> apply_masking(const vector<complex<double>>& X, double masking_factor=0.15){
    vector<complex<double>> Y = X;
    double max_mag = 0.0;
    for(auto c : X) max_mag = max(max_mag, abs(c));
    for(auto &c : Y){
        if(abs(c) < masking_factor*max_mag) c = complex<double>(0,0);
    }
    return Y;
}

/* ===== Lossless compression: Run-length encoding (RLE) ===== */
vector<pair<double,int>> lossless_rle(const vector<double>& x){
    vector<pair<double,int>> rle;
    if(x.empty()) return rle;
    double val = x[0];
    int count = 1;
    for(size_t i=1;i<x.size();i++){
        if(x[i]==val) count++;
        else {
            rle.push_back({val,count});
            val = x[i];
            count = 1;
        }
    }
    rle.push_back({val,count});
    return rle;
}

/* ===== Lossy compression: uniform quantization ===== */
vector<double> lossy_quantize(const vector<double>& x, int levels){
    vector<double> y(x.size());
    double min_val = *min_element(x.begin(), x.end());
    double max_val = *max_element(x.begin(), x.end());
    double q_step = (max_val-min_val)/(levels-1);
    for(size_t i=0;i<x.size();i++){
        y[i] = round((x[i]-min_val)/q_step)*q_step + min_val;
    }
    return y;
}

int main(){
    cout << "=== Compression Demo ===" << endl;

    // Example signal (numeric or audio-like)
    vector<double> signal(64);
    for(size_t i=0;i<signal.size();i++){
        signal[i] = sin(2*PI*5*i/signal.size()) + 0.5*sin(2*PI*20*i/signal.size());
    }
    export_csv("original_signal.csv", signal);

    // ===== Lossless compression demonstration =====
    auto rle = lossless_rle(signal);
    cout << "\nLossless RLE compression result (first 10 pairs):" << endl;
    for(size_t i=0;i<min(size_t(10), rle.size());i++){
        cout << "(" << rle[i].first << "," << rle[i].second << ") ";
    }
    cout << endl;

    // ===== Lossy compression demonstration =====
    auto quantized = lossy_quantize(signal, 8);
    export_csv("quantized_signal.csv", quantized);
    cout << "\nLossy quantization with 8 levels done." << endl;

    // ===== Psychoacoustic-inspired masking =====
    auto X = DFT(signal);
    auto X_masked = apply_masking(X, 0.15);
    auto reconstructed = IDFT(X_masked);
    export_csv("masked_signal.csv", reconstructed);
    cout << "\nPsychoacoustic masking completed." << endl;

    cout << "\nCSV files exported: original, quantized, masked/reconstructed." << endl;
    return 0;
}

/*
# Python plotting code

import pandas as pd
import matplotlib.pyplot as plt

orig = pd.read_csv('original_signal.csv')['value']
quant = pd.read_csv('quantized_signal.csv')['value']
masked = pd.read_csv('masked_signal.csv')['value']

plt.figure(figsize=(12,6))
plt.plot(orig, label='Original Signal')
plt.plot(quant, label='Lossy Quantized', linestyle='--')
plt.plot(masked, label='Masked/Reconstructed', linestyle=':')
plt.title("Lossless vs Lossy Compression and Psychoacoustic Masking")
plt.xlabel("Sample Index")
plt.ylabel("Amplitude")
plt.grid(True)
plt.legend()
plt.show()
*/
