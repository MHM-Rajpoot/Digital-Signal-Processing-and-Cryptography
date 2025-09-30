#include "armadillo-15.0.2\include\armadillo"
#include "sigpack-1.2.7\sigpack\sigpack.h"
#include <iostream>

using namespace sp;
using namespace arma;
using namespace std;

int main() {

    vec x = linspace(0, 2*M_PI, 10);
    vec y = sin(x);
    cout << "Relative include test: " << y.t();

    return 0;
}
