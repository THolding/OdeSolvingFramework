#pragma once
#include <iostream>
#include <cmath>
#include "utilities.hpp"
#include "model_definition.hpp"

class SIRSModel : public ModelDefinition
{
private:
    //
public:
    SIRSModel()
    {  }


    //Parameters: ...
    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        const double beta = params[0];
        const double sigma = params[1];
        const double mu = params[2];
        const double alpha = params[3];
        const double N = 1.0;

        const double S = currentValues[0];
        const double I = currentValues[1];
        const double R = currentValues[2];

        ////Calculate dxdt for each compartment.
        output[0] = -beta*S*I +alpha*R +mu*N -mu*S;  //Susceptible compartment.
        output[1] =  beta*S*I -sigma*I -mu*I; //Infected.
        output[2] =  sigma*I  -alpha*R -mu*R; //Recovered.
    }
};
