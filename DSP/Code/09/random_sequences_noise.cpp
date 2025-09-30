/*
 * random_sequences_noise.cpp
 *
 * Demonstration of random sequences and noise:
 *  - Random variables
 *  - Stationary processes
 *  - Autocorrelation and cross-correlation
 *  - Deterministic crosscorrelation sequences
 *  - Filtered random sequences
 *  - White noise generation
 *  - Exponential averaging
 *
 * CSV export for MATLAB/Python/GNUplot.
 *
 * Inspired by Dr Markus Kuhn, Computer Laboratory, University of Cambridge.
 *
 * Compile with: g++ -std=c++11 -o random_sequences_noise random_sequences_noise.cpp -lm
 */

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
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

/* ===== Random sequences ===== */
vector<double> generate_uniform(int N, double a=0.0, double b=1.0) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(a,b);
    vector<double> x(N);
    for(int i=0;i<N;i++) x[i] = dis(gen);
    return x;
}

vector<double> generate_gaussian(int N, double mu=0.0, double sigma=1.0) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dis(mu,sigma);
    vector<double> x(N);
    for(int i=0;i<N;i++) x[i] = dis(gen);
    return x;
}

/* ===== Autocorrelation ===== */
vector<double> autocorrelation(const vector<double>& x) {
    int N = x.size();
    vector<double> r(N,0.0);
    for(int k=0;k<N;k++){
        for(int n=0;n<N-k;n++){
            r[k] += x[n]*x[n+k];
        }
        r[k] /= N;
    }
    return r;
}

/* ===== Cross-correlation ===== */
vector<double> crosscorrelation(const vector<double>& x, const vector<double>& y) {
    int N = x.size();
    vector<double> r(N,0.0);
    for(int k=0;k<N;k++){
        for(int n=0;n<N-k;n++){
            r[k] += x[n]*y[n+k];
        }
        r[k] /= N;
    }
    return r;
}

/* ===== Deterministic cross-correlation sequence ===== */
vector<double> deterministic_sequence(int N) {
    vector<double> seq(N);
    for(int n=0;n<N;n++)
        seq[n] = (n%2==0?1.0:-1.0); // simple alternating +1/-1 sequence
    return seq;
}

/* ===== Filtered random sequence ===== */
vector<double> lowpass_filter(const vector<double>& x, double alpha) {
    vector<double> y(x.size(),0.0);
    y[0] = x[0];
    for(size_t n=1;n<x.size();n++)
        y[n] = alpha*x[n] + (1-alpha)*y[n-1];
    return y;
}

/* ===== Exponential averaging ===== */
vector<double> exponential_average(const vector<double>& x, double alpha) {
    vector<double> y(x.size(),0.0);
    y[0] = x[0];
    for(size_t n=1;n<x.size();n++)
        y[n] = alpha*x[n] + (1-alpha)*y[n-1];
    return y;
}

/* ===== Main demonstration ===== */
int main(){
    cout<<"=== Random Sequences and Noise Demonstration ==="<<endl;

    int N = 128;

    // Random variables: uniform and Gaussian
    auto x_uniform = generate_uniform(N);
    auto x_gaussian = generate_gaussian(N);

    export_csv("uniform_sequence.csv",x_uniform);
    export_csv("gaussian_sequence.csv",x_gaussian);

    // Stationary process: simple filtered Gaussian
    auto x_filtered = lowpass_filter(x_gaussian,0.1);
    export_csv("filtered_gaussian.csv",x_filtered);

    // Autocorrelation
    auto r_auto = autocorrelation(x_gaussian);
    export_csv("autocorrelation.csv",r_auto);

    // Cross-correlation with deterministic sequence
    auto det_seq = deterministic_sequence(N);
    auto r_cross = crosscorrelation(x_gaussian,det_seq);
    export_csv("crosscorrelation.csv",r_cross);

    // Exponential averaging
    auto x_expavg = exponential_average(x_gaussian,0.1);
    export_csv("exponential_average.csv",x_expavg);

    cout<<"All CSV files exported. Plot using Python/MATLAB/GNUplot."<<endl;

    return 0;
}

/*
# Python plotting code (save as plot_random_sequences.py)
import pandas as pd
import matplotlib.pyplot as plt

def plot_sequence(file,title):
    data = pd.read_csv(file)
    plt.figure(figsize=(6,3))
    plt.plot(data['n'], data.iloc[:,1])
    plt.title(title)
    plt.xlabel('n')
    plt.ylabel('Amplitude')
    plt.grid(True)

plot_sequence('uniform_sequence.csv','Uniform Random Sequence')
plot_sequence('gaussian_sequence.csv','Gaussian Random Sequence')
plot_sequence('filtered_gaussian.csv','Filtered Gaussian Sequence')
plot_sequence('autocorrelation.csv','Autocorrelation')
plot_sequence('crosscorrelation.csv','Cross-correlation with Deterministic Sequence')
plot_sequence('exponential_average.csv','Exponential Average of Gaussian Sequence')

plt.show()
*/
