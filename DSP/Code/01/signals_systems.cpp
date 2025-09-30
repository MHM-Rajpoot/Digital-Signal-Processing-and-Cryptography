/* 
 * signals_systems.cpp
 *
 * Demonstration of discrete sequences and systems, their types
 * and properties, including linear time-invariant (LTI) systems,
 * convolution, linearity, and time-invariance.
 *
 * Includes CSV export for plotting in MATLAB, Python, or GNUplot.
 *
 * Written in the style of Computer Laboratory examples
 * (inspired by Dr Markus Kuhn's teaching material).
 *
 * Compile with: g++ -std=c++11 -o signals_systems signals_systems.cpp
 */

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>   // for CSV export

using namespace std;

/* Utility: print a discrete sequence */
void print_sequence(const vector<double> &x, const string &name) {
    cout << name << " = { ";
    for (size_t n = 0; n < x.size(); ++n) {
        cout << fixed << setprecision(2) << x[n];
        if (n < x.size() - 1) cout << ", ";
    }
    cout << " }" << endl;
}

/* Utility: export a sequence to CSV */
void export_to_csv(const vector<double> &x, const string &filename, const string &colname="x[n]") {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: could not open file " << filename << endl;
        return;
    }

    file << "n," << colname << "\n";
    for (size_t n = 0; n < x.size(); ++n) {
        file << n << "," << x[n] << "\n";
    }

    file.close();
    cout << "Exported sequence to " << filename << endl;
}

/* Example 1: Define basic discrete sequences */
vector<double> unit_impulse(int N, int k) {
    // δ[n-k] of length N
    vector<double> x(N, 0.0);
    if (k >= 0 && k < N) x[k] = 1.0;
    return x;
}

vector<double> unit_step(int N, int k) {
    // u[n-k] of length N
    vector<double> x(N, 0.0);
    for (int n = k; n < N; ++n) {
        x[n] = 1.0;
    }
    return x;
}

/* Example 2: Convolution implementation */
vector<double> convolution(const vector<double> &x, const vector<double> &h) {
    int N = x.size();
    int M = h.size();
    vector<double> y(N + M - 1, 0.0);

    for (int n = 0; n < N + M - 1; ++n) {
        for (int k = 0; k < N; ++k) {
            if (n - k >= 0 && n - k < M) {
                y[n] += x[k] * h[n - k];
            }
        }
    }
    return y;
}

/* Example 3a: Simple LTI system (moving average filter) */
vector<double> moving_average(const vector<double> &x, int L) {
    int N = x.size();
    vector<double> y(N, 0.0);

    for (int n = 0; n < N; ++n) {
        double sum = 0.0;
        for (int k = 0; k < L; ++k) {
            if (n - k >= 0) sum += x[n - k];
        }
        y[n] = sum / L;
    }
    return y;
}

/* Example 3b: Non-time-invariant system (y[n] = n * x[n]) */
vector<double> time_varying_system(const vector<double> &x, int dummy) {
    int N = x.size();
    vector<double> y(N, 0.0);

    for (int n = 0; n < N; ++n) {
        y[n] = n * x[n]; // depends explicitly on time index
    }
    return y;
}

/* Example 4: Check linearity of a system */
bool check_linearity(const vector<double> &x1, const vector<double> &x2,
                     double a, double b, 
                     vector<double> (*system)(const vector<double>&, int), int L) {
    int N = x1.size();
    vector<double> combined(N);
    for (int i = 0; i < N; ++i)
        combined[i] = a * x1[i] + b * x2[i];

    vector<double> lhs = system(combined, L);

    vector<double> s1 = system(x1, L);
    vector<double> s2 = system(x2, L);
    vector<double> rhs(N);
    for (int i = 0; i < N; ++i)
        rhs[i] = a * s1[i] + b * s2[i];

    for (int i = 0; i < N; ++i)
        if (abs(lhs[i] - rhs[i]) > 1e-6) return false;

    return true;
}

/* Example 5: Check time-invariance of a system */
bool check_time_invariance(const vector<double> &x, int shift,
                           vector<double> (*system)(const vector<double>&, int), int L) {
    int N = x.size();
    
    // Shift input
    vector<double> x_shifted(N, 0.0);
    for (int n = shift; n < N; ++n)
        x_shifted[n] = x[n - shift];

    // Apply system
    vector<double> y_original = system(x, L);
    vector<double> y_shifted_input = system(x_shifted, L);

    // Shift output
    vector<double> y_shifted_output(N, 0.0);
    for (int n = shift; n < N; ++n)
        y_shifted_output[n] = y_original[n - shift];

    // Compare both
    for (int i = 0; i < N; ++i)
        if (abs(y_shifted_input[i] - y_shifted_output[i]) > 1e-6) return false;

    return true;
}

int main() {
    cout << "=== Signals and Systems Demonstration ===\n" << endl;

    // Basic sequences
    auto delta = unit_impulse(8, 0);
    auto step  = unit_step(8, 2);

    print_sequence(delta, "Unit impulse δ[n]");
    print_sequence(step, "Unit step u[n-2]");

    export_to_csv(delta, "delta.csv", "δ[n]");
    export_to_csv(step, "step.csv", "u[n-2]");

    // Convolution example
    vector<double> x = {1.0, 2.0, 3.0};
    vector<double> h = {0.5, 0.5};  // Simple averaging filter
    auto y = convolution(x, h);

    print_sequence(x, "x[n]");
    print_sequence(h, "h[n]");
    print_sequence(y, "y[n] = x[n] * h[n]");

    export_to_csv(x, "x.csv", "x[n]");
    export_to_csv(h, "h.csv", "h[n]");
    export_to_csv(y, "y_conv.csv", "y[n]");

    // Moving average as an LTI system
    vector<double> input = {1, 2, 0, -1, 3, 2};
    auto filtered = moving_average(input, 3);
    print_sequence(input, "Input");
    print_sequence(filtered, "Moving average output (LTI system)");

    export_to_csv(input, "input.csv", "input[n]");
    export_to_csv(filtered, "moving_avg.csv", "y[n] (moving average)");

    // Linearity check
    vector<double> x1 = {1, 2, 3};
    vector<double> x2 = {0, -1, 1};
    bool is_linear = check_linearity(x1, x2, 2.0, -1.5, moving_average, 3);
    cout << "\nLinearity of moving average system: "
         << (is_linear ? "YES" : "NO") << endl;

    // Time-invariance check: moving average
    bool is_time_invariant = check_time_invariance(input, 2, moving_average, 3);
    cout << "Time-invariance of moving average system: "
         << (is_time_invariant ? "YES" : "NO") << endl;

    // Time-invariance check: time-varying system y[n] = n*x[n]
    bool tv_invariant = check_time_invariance(input, 2, time_varying_system, 0);
    cout << "Time-invariance of y[n] = n*x[n] system: "
         << (tv_invariant ? "YES" : "NO") << endl;

    return 0;
}

/* =====================================================================
   Python plotting code (save separately as signals_systems.cpp)
   =====================================================================

import pandas as pd
import matplotlib.pyplot as plt

def plot_sequence(filename, colname, title):
    data = pd.read_csv(filename)
    plt.figure(figsize=(6,4))
    plt.stem(data['n'], data[colname], basefmt='k-', use_line_collection=True)
    plt.title(title)
    plt.xlabel('n')
    plt.ylabel(colname)
    plt.grid(True)

# Plot unit impulse δ[n]
plot_sequence('delta.csv', 'δ[n]', 'Unit Impulse δ[n]')

# Plot unit step u[n-2]
plot_sequence('step.csv', 'u[n-2]', 'Unit Step u[n-2]')

# Plot x[n], h[n], and y[n] from convolution
plot_sequence('x.csv', 'x[n]', 'Signal x[n]')
plot_sequence('h.csv', 'h[n]', 'Impulse Response h[n]')
plot_sequence('y_conv.csv', 'y[n]', 'Convolution y[n] = x[n] * h[n]')

# Plot moving average system input/output
plot_sequence('input.csv', 'input[n]', 'Input Signal to Moving Average')
plot_sequence('moving_avg.csv', 'y[n] (moving average)', 'Moving Average Output (LTI System)')

plt.tight_layout()
plt.show()
*/
