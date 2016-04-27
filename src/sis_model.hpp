#pragma once
#include <iostream>
#include <cmath>
#include "utilities.hpp"
#include "model_definition.hpp"

class SISModel : public ModelDefinition
{
private:
    //
public:
    SISModel()
    {  }


    //Parameters: ...
    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        const double beta = params[0];
        const double sigma = params[1];
        const double mu = params[2];
        const double N = 1.0;

        const double S = currentValues[0];
        const double I = currentValues[1];


        ////Calculate dxdt for each compartment.
        output[0] = -beta*S*I +sigma*I +mu*N -mu*S;  //Susceptible compartment.
        output[1] =  beta*S*I -sigma*I -mu*I; //Infected.

    }
};
