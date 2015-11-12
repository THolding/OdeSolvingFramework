#pragma once
#include <vector>
#include <stdio.h> //temp, replace with C++ style file io

class ModelDriver;

//TODO: Include seasonality as a subclass.
//Class encapsulating the Runge-Kutta Cash-Karp method with an adaptive time step for numerically solving differential equations.
class RungeKuttaCKSolver
{
private:
    //Cash-Karp parameters for embedded Runge-Kutta:
    static constexpr double b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42=-0.9, b43=1.2, b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0,
        b54=35.0/27.0, b61=1631.0/55296, b62=175.0/512.0, b63=575.0/13824.0, b64=44275.0/110592, b65=253.0/4096.0;
    static constexpr double c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0;
    static constexpr double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, dc5=-277.00/14336, dc6=c6-0.25;

    //Default values:
    static constexpr double defaultPrecision = 1E-10; //UNUSED because we have non-adaptive timestep.

    //Solver specific fields:
    ModelDriver* model;
    double precision; //Required precision. Used to guide adaptive timestep. UNUSED because we have non-adaptive timestep.

public:
    RungeKuttaCKSolver(ModelDriver* _model);
    //double[] getDerivatives(double* _yCurrent); //Gets the derivatives for each equation
    void calc(const std::vector<double>& _yCurrent, const std::vector<double>& _derivs, const double _timestep, std::vector<double>& _yOut, std::vector<double>& _yErr); //Performs a single 'dumb' trial step
    void step(std::vector<double>& _yCurrent, const std::vector<double> &_derivs, double *_stepsize, double *_nextStepsize/*, const std::vector<double> &_yScale*/); //Stepper (manages adaptive time steps)
    void solve(); //Driver
};
