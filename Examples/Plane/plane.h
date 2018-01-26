#ifndef PLANE_H
#define PLANE_H

#include <cmath>
#include <extended.h>

class PlaneFilter : public Kalman::Extended<double>
{
public:
    PlaneFilter();

protected:
    void makeA() override;
    void makeH() override;
    void makeV() override;
    void makeR() override;
    void makeW() override;
    void makeQ() override;
    void makeProcess() override;
    void makeMeasure() override;

protected:
    double Period;
    double Mass;
    double Drag;
    double Lift;
    double Gravity;
};

#endif
