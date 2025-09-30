/*
 * digital_modulation_demo.cpp
 *
 * Demonstration of digital modulation techniques:
 *  - IQ representation of AM, FM, MSK, QAM, OFDM
 *  - Symbol detection
 *  - Clock recovery (simple sampling)
 *  - Matched filtering
 *
 * CSV export for MATLAB/Python/GNUplot.
 *
 * Inspired by Dr Markus Kuhn's Computer Laboratory material.
 *
 * Compile with: g++ -std=c++11 -o digital_modulation_demo digital_modulation_demo.cpp -lm
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>

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

/* ===== Generate random bit sequence ===== */
vector<int> random_bits(int N) {
    vector<int> bits(N);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0,1);
    for(int i=0;i<N;i++) bits[i] = dis(gen);
    return bits;
}

/* ===== IQ Representation and Modulation ===== */
vector<complex<double>> am_modulate(const vector<double>& m, double A, double fc, double fs) {
    vector<complex<double>> s(m.size());
    double Ts = 1.0/fs;
    for(size_t n=0;n<m.size();n++)
        s[n] = A*(1+m[n])*exp(complex<double>(0,2*M_PI*fc*n*Ts));
    return s;
}

vector<complex<double>> fm_modulate(const vector<double>& m, double A, double fc, double kf, double fs) {
    vector<complex<double>> s(m.size());
    double Ts = 1.0/fs;
    double phase = 0.0;
    for(size_t n=0;n<m.size();n++) {
        phase += 2*M_PI*kf*m[n]*Ts;
        s[n] = A*exp(complex<double>(0,2*M_PI*fc*n*Ts + phase));
    }
    return s;
}

vector<complex<double>> qam_modulate(const vector<int>& bits, int M, double fc, double fs) {
    vector<complex<double>> s(bits.size()/2);
    double Ts = 1.0/fs;
    int idx=0;
    for(size_t i=0;i<bits.size();i+=2) {
        double I = (bits[i]==0?-1:1);
        double Q = (bits[i+1]==0?-1:1);
        s[idx++] = complex<double>(I,Q)*exp(complex<double>(0,2*M_PI*fc*(idx-1)*Ts));
    }
    return s;
}

/* Simple rectangular pulse shaping for MSK (1 bit per symbol) */
vector<complex<double>> msk_modulate(const vector<int>& bits, double fc, double fs) {
    vector<complex<double>> s(bits.size());
    double Ts = 1.0/fs;
    double phase = 0.0;
    for(size_t n=0;n<bits.size();n++) {
        phase += M_PI*bits[n]; // 0->0, 1->pi phase jump
        s[n] = exp(complex<double>(0,2*M_PI*fc*n*Ts + phase));
    }
    return s;
}

/* OFDM generation: simple 4 subcarriers, no cyclic prefix */
vector<complex<double>> ofdm_modulate(const vector<int>& bits, int N_subcarriers, double fc, double fs) {
    int N_sym = bits.size()/N_subcarriers;
    vector<complex<double>> s(N_sym*N_subcarriers,0);
    double Ts = 1.0/fs;
    for(int k=0;k<N_sym;k++)
        for(int n=0;n<N_subcarriers;n++){
            double I = (bits[k*N_subcarriers+n]==0?-1:1);
            s[k*N_subcarriers+n] = I*exp(complex<double>(0,2*M_PI*(fc+n*0.01)*((k*N_subcarriers+n)*Ts)));
        }
    return s;
}

/* Simple matched filter (correlator) */
vector<double> matched_filter(const vector<complex<double>>& s, const vector<complex<double>>& template_signal){
    int N = s.size()-template_signal.size()+1;
    vector<double> output(N,0.0);
    for(int n=0;n<N;n++){
        complex<double> sum=0;
        for(size_t k=0;k<template_signal.size();k++)
            sum += s[n+k]*conj(template_signal[k]);
        output[n] = abs(sum);
    }
    return output;
}

/* Main demonstration */
int main(){
    cout<<"=== Digital Modulation Demonstration ==="<<endl;

    int N_bits = 64;
    double fs = 1000.0;
    double fc = 100.0;

    // Random bits
    auto bits = random_bits(N_bits);

    // Map bits to +/-1 for AM/FM/MSK
    vector<double> bits_double(bits.begin(), bits.end());
    for(auto &v: bits_double) v = 2*v -1;

    // AM
    auto s_am = am_modulate(bits_double,1.0,fc,fs);
    export_csv("am_signal.csv",s_am);

    // FM
    auto s_fm = fm_modulate(bits_double,1.0,fc,0.1,fs);
    export_csv("fm_signal.csv",s_fm);

    // MSK
    auto s_msk = msk_modulate(bits,fc,fs);
    export_csv("msk_signal.csv",s_msk);

    // QAM (4-QAM)
    auto s_qam = qam_modulate(bits,4,fc,fs);
    export_csv("qam_signal.csv",s_qam);

    // OFDM
    auto s_ofdm = ofdm_modulate(bits,4,fc,fs);
    export_csv("ofdm_signal.csv",s_ofdm);

    // Matched filter example: detect first 8-bit sequence in AM signal
    vector<complex<double>> template_signal(s_am.begin(), s_am.begin()+8);
    auto mf_output = matched_filter(s_am,template_signal);
    export_csv("matched_filter_output.csv",mf_output);

    cout<<"All CSV files exported. Plot in Python/MATLAB/GNUplot."<<endl;

    return 0;
}

/*
# Python plotting code for digital modulation (save as digital_mod_plot.py)
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
plt.stem(mf['n'], mf['value'], use_line_collection=True)
plt.title('Matched Filter Output')
plt.xlabel('n')
plt.ylabel('Correlation Magnitude')
plt.grid(True)

plt.show()
*/
