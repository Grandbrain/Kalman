#include <iostream>
#include <fstream>
#include <iomanip>
#include "plane.h"

using namespace std;
using namespace Kalman;

int main()
{
    const unsigned count = 500;
    const unsigned n = 4; // states
    const unsigned m = 2; // measures

    const double raw[] = { 100.0 * 100.0, 0.0, 0.0, 0.0,
                           0.0, 10.0 * 10.0, 0.0, 0.0,
                           0.0, 0.0, 25.0 * 25.0, 0.0,
                           0.0, 0.0, 0.0, 10.0 * 10.0 };

    ifstream input;
    PlaneFilter filter;

    Vector<double> x(n);
    Vector<double> z(m);
    Matrix<double> p(n, n, raw);

    // read the inputs vector and the measures matrix
    input.open("data.dat", ifstream::in);

    double fs[count] = {0.0};
    double ms[count * 2] = {0.0};

    for (double& fi : fs) input >> fi;
    for (double& mi : ms) input >> mi;

    Vector<double> f(count, fs);
    Matrix<double> measure(m, count, ms);

    // initialization
    x(0) = cos(measure(0,0)) * measure(1,0);
    x(1) = 60;
    x(2) = sin(measure(0,0)) * measure(1,0);
    x(3) = 0;

    filter.init(x, p);

    // process
    for (unsigned i = 1; i < count; ++i)
    {
        for(unsigned j = 0; j < m; ++j)
            z(j) = measure(j,i);

        Vector<double> u(1, f(i));

        filter.step(u, z);

        Vector<double> result = filter.getX();

        cout << result(0) << ' ' << result(1) << ' '
             << result(2) << ' ' << result(3) << endl;
    }

    return 0;
}
