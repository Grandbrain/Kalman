#include "plane.h"

PlaneFilter::PlaneFilter()
{
    setSizeX(4);
    setSizeU(1);
    setSizeW(2);
    setSizeZ(2);
    setSizeV(2);

    Period = 0.2;
    Mass = 1000;
    Drag = 0.35;
    Lift = 3.92;
    Gravity = 9.8;
}

void PlaneFilter::makeA()
{
    A(0,0) = 1.0;
    A(0,1) = Period - Period * Period * Drag / Mass * x(1);
    A(0,2) = 0.0;
    A(0,3) = 0.0;

    A(1,0) = 0.0;
    A(1,1) = 1 - 2 * Period * Drag / Mass * x(1);
    A(1,2) = 0.0;
    A(1,3) = 0.0;

    A(2,0) = 0.0;
    A(2,1) = Period * Period * Lift / Mass * x(1);
    A(2,2) = 1.0;
    A(2,3) = Period;

    A(3,0) = 0.0;
    A(3,1) = 2 * Period * Lift / Mass * x(1);
    A(3,2) = 0.0;
    A(3,3) = 1.0;
}

void PlaneFilter::makeW()
{
    W(0,0) = 0.0;
    W(0,1) = 0.0;

    W(1,0) = 1.0;
    W(1,1) = 0.0;

    W(2,0) = 0.0;
    W(2,1) = 0.0;

    W(3,0) = 0.0;
    W(3,1) = 1.0;
}

void PlaneFilter::makeQ()
{
    Q(0,0) = 0.01 * 0.01;
    Q(0,1) = 0.01 * 0.01/10.0;

    Q(1,0) = 0.01 * 0.01 / 10.0;
    Q(1,1) = 0.01 * 0.01;
}

void PlaneFilter::makeH()
{
    H(0,0) = -x(2) / (x(0) * x(0) + x(2) * x(2));
    H(0,1) = 0.0;
    H(0,2) = x(0) / (x(0) * x(0) + x(2) * x(2));
    H(0,3) = 0.0;

    H(1,0) = x(0) / sqrt(x(0) * x(0) + x(2) * x(2));
    H(1,1) = 0.0;
    H(1,2) = x(2) / sqrt(x(0) * x(0) + x(2) * x(2));
    H(1,3) = 0.0;
}

void PlaneFilter::makeV()
{
    V(0,0) = 1.0;
    V(0,1) = 0.0;

    V(1,0) = 0.0;
    V(1,1) = 1.0;
}

void PlaneFilter::makeR()
{
    R(0,0) = 0.01 * 0.01;
    R(0,1) = 0.0;

    R(1,0) = 0.0;
    R(1,1) = 50 * 50;
}

void PlaneFilter::makeProcess()
{
    Kalman::Vector<double> x_(x.Size());
    x_(0) = x(0) + x(1) * Period + (Period * Period) / 2 * (u(0) / Mass -
            Drag / Mass * x(1) * x(1));
    x_(1) = x(1) + (u(0) / Mass - Drag / Mass * x(1) * x(1)) * Period;
    x_(2) = x(2) + x(3) * Period + (Period * Period) / 2 *
                                   (Lift / Mass * x(1) * x(1) - Gravity);
    x_(3) = x(3) + (Lift / Mass * x(1) * x(1) - Gravity) * Period;
    x.Swap(x_);
}

void PlaneFilter::makeMeasure()
{
    z(0) = atan2(x(2), x(0));
    z(1) = sqrt(x(0) * x(0) + x(2) * x(2));
}
