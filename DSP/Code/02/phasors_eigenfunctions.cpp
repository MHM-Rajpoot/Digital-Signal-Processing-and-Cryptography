/* 
 * phasors_eigenfunctions.cpp
 *
 * Demonstration of phasors, eigenfunctions of LTI systems,
 * review of complex arithmetic, and applications in 
 * electronics, optics, and acoustics.
 *
 * Includes CSV export for plotting in MATLAB, Python, or GNUplot.
 *
 * Written in the style of Computer Laboratory examples
 * (inspired by Dr Markus Kuhn's teaching material).
 *
 * Compile with: g++ -std=c++11 -o phasors_eigenfunctions phasors_eigenfunctions.cpp
 */

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

/* Utility: print a complex number in polar form */
void print_polar(const complex<double>& z, const string &name) {
    double mag = abs(z);
    double phase = arg(z); // radians
    cout << name << " = " << fixed << setprecision(3)
         << mag << " ∠ " << phase << " rad" << endl;
}

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

/* Example 1: Review of complex arithmetic */
void complex_arithmetic_demo() {
    cout << "\n=== Complex Arithmetic Review ===" << endl;

    complex<double> a(3, 4);  // 3 + j4
    complex<double> b(1, -2); // 1 - j2

    cout << "a = " << a << ", b = " << b << endl;
    cout << "a + b = " << a + b << endl;
    cout << "a * b = " << a * b << endl;
    cout << "conj(a) = " << conj(a) << endl;
    print_polar(a, "|a| ∠ arg(a)");
}

/* Example 2: Phasors */
void phasor_demo() {
    cout << "\n=== Phasors ===" << endl;

    double A = 2.0;      // amplitude
    double omega = M_PI; // angular frequency = pi rad/s
    double phi = M_PI/4; // phase shift = 45 deg

    cout << "Signal: x(t) = " << A << " cos(ωt + φ)" << endl;

    // Represent as phasor
    complex<double> phasor = polar(A, phi); // A∠φ
    print_polar(phasor, "Phasor representation");

    // Example: multiply by system gain H = 0.5∠-30°
    complex<double> H = polar(0.5, -M_PI/6);
    complex<double> Y = H * phasor;
    print_polar(H, "System response H");
    print_polar(Y, "Output phasor Y");

    // Export phasor data
    vector<complex<double>> data = {phasor, H, Y};
    export_csv("phasors.csv", data);
}

/* Example 3: Eigenfunctions of LTI systems */
void eigenfunction_demo() {
    cout << "\n=== Eigenfunctions of LTI Systems ===" << endl;

    // Complex exponential input x[n] = e^(jωn)
    double omega = M_PI / 3; // frequency
    vector<complex<double>> x;
    for (int n = 0; n < 6; ++n) {
        x.push_back(exp(complex<double>(0, omega * n)));
    }

    cout << "Input x[n] = e^(jωn): ";
    for (auto &val : x) cout << val << " ";
    cout << endl;

    // Example LTI system: FIR filter h[n] = {1, 0.5}
    vector<double> h = {1.0, 0.5};
    vector<complex<double>> y(x.size() + h.size() - 1, 0);

    for (int n = 0; n < (int)y.size(); ++n) {
        for (int k = 0; k < (int)x.size(); ++k) {
            if (n - k >= 0 && n - k < (int)h.size())
                y[n] += x[k] * h[n - k];
        }
    }

    cout << "Output y[n] = convolution result: ";
    for (auto &val : y) cout << val << " ";
    cout << endl;

    // Theoretical eigenvalue H(e^jω)
    complex<double> H_eval = 0.0;
    for (int n = 0; n < (int)h.size(); ++n)
        H_eval += h[n] * exp(complex<double>(0, -omega * n));

    print_polar(H_eval, "Eigenvalue H(e^jω)");

    // Export signals
    export_csv("eigen_x.csv", x);
    export_csv("eigen_y.csv", y);
    vector<complex<double>> H_data = {H_eval};
    export_csv("eigen_H.csv", H_data);
}

/* Example 4: Applications */
void applications_demo() {
    cout << "\n=== Applications ===" << endl;

    // Electronics: AC circuit (R + jX)
    complex<double> Z_elec(10, 5); // impedance = 10 + j5 Ω
    complex<double> I(2, 0);       // current phasor 2∠0
    complex<double> V = Z_elec * I;
    print_polar(V, "Electronics: Voltage phasor across R-L circuit");

    // Optics: Plane wave E = E0 exp(j(kz - ωt))
    double k = 2 * M_PI / 500e-9;  // wavenumber for λ=500nm
    complex<double> E0(1.0, 0.0);  // unit amplitude
    complex<double> E = E0 * exp(complex<double>(0, k * 1e-6));
    print_polar(E, "Optics: Electric field phasor at z=1 µm");

    // Acoustics: Sound pressure phasor
    complex<double> P0 = polar(1.0, M_PI/6); // 1 Pa ∠30°
    complex<double> Z_ac(400, 0);            // acoustic impedance (real)
    complex<double> U = P0 / Z_ac;           // volume velocity
    print_polar(U, "Acoustics: Volume velocity phasor");

    // Export application results
    vector<complex<double>> app_data = {V, E, U};
    export_csv("applications.csv", app_data);
}

int main() {
    cout << "=== Phasors and Eigenfunctions Demonstration ===\n" << endl;

    complex_arithmetic_demo();
    phasor_demo();
    eigenfunction_demo();
    applications_demo();
    
    return 0;
}

/* =====================================================================
   Python plotting code (save separately as plot_phasors_eigen.py)
   =====================================================================

import pandas as pd
import matplotlib.pyplot as plt

def plot_time_signal(filename, title):
    data = pd.read_csv(filename)
    plt.figure(figsize=(6,4))
    plt.stem(data['n'], data['real'], linefmt='b-', markerfmt='bo', basefmt='k-')
    plt.title(title + " (Real Part)")
    plt.xlabel("n")
    plt.ylabel("Value")
    plt.grid(True)

def plot_phasors(filename, title):
    data = pd.read_csv(filename)
    plt.figure(figsize=(6,6))
    plt.quiver([0]*len(data), [0]*len(data),
               data['real'], data['imag'],
               angles='xy', scale_units='xy', scale=1, color=['r','g','b'])
    for i, txt in enumerate(data.index):
        plt.text(data['real'][i]*1.05, data['imag'][i]*1.05, f"{i}")
    plt.title(title + " (Phasors)")
    plt.xlabel("Re")
    plt.ylabel("Im")
    plt.grid(True)
    plt.axis('equal')

# Plotting
plot_phasors("phasors.csv", "Phasor Demo")
plot_time_signal("eigen_x.csv", "Eigenfunction Input x[n]")
plot_time_signal("eigen_y.csv", "Eigenfunction Output y[n]")
plot_phasors("applications.csv", "Applications (Electronics, Optics, Acoustics)")

plt.show()

===================================================================== */
