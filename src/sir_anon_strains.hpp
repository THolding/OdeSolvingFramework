#pragma once
#include <iostream>
#include <cmath>
#include "utilities.hpp"
#include "model_definition.hpp"

class SIRAnonStrains : public ModelDefinition
{
private:
    const unsigned int numStrains;
public:
    SIRAnonStrains(const unsigned int _numStrains) : numStrains(_numStrains)
    {  }
    //Parameters: numStrains+1 betas (final beta is typically 0, indicating total immunity).
    //            numStrains sigmas.
    //            mu at params.end()-1.

    // Compartments as follows: S0 -> Sn-1, I0 -> In-1, R.
    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        //Define parameter constants for easier readability.
        const std::vector<double> betas(params.begin(), params.begin()+numStrains);
        const std::vector<double> sigmas(params.begin()+numStrains, params.begin()+(2*numStrains));
        const double mu = params.back();
        const double N = 1.0;
        double totalInfected = 0;
        for (unsigned int i=0; i<numStrains; i++)
            totalInfected += currentValues[numStrains+i];
        const double R = currentValues[numStrains*2];


        //Calculate dxdt for the first susceptible compartment.
        output[0] = -betas[0]*currentValues[0]*totalInfected +mu*N -mu*currentValues[0];

        //Other susceptible compartments
        for (unsigned int i=1; i<numStrains; i++)
        {
            const double I_im1 = currentValues[numStrains+i-1];
            output[i] = -betas[i]*currentValues[i]*totalInfected +sigmas[i-1]*I_im1 -mu*currentValues[i];
        }

        //Infected compartments
        for (unsigned int i=0; i<numStrains; i++)
        {
            const double I_i = currentValues[numStrains+i];
            output[numStrains+i] = betas[i]*currentValues[i]*totalInfected -sigmas[i]*I_i -mu*I_i;
        }

        //Recovered compartment
        output[numStrains*2] = +sigmas[numStrains-1]*currentValues[(numStrains*2)-1] -mu*R;

        /*
        //Calculate dxdt for the first susceptible compartment.
        output[0] = -betas[0]*currentValues[0]*totalInfected +mu*N -mu*currentValues[0]; //First compartment is a special case because it has no inflow flux from an infected compartment.

        //Calculate dxdt for middle susceptible compartments.
        for (unsigned int i=1; i<numStrains; i++)
        {
            const double S_i = currentValues[i];
            const double I_im1 = currentValues[numStrains+i-1];

            //Note that the first 50% of compartments are taken to be the susceptible compartments, while the second 50% are infected, with the recovered compartment in the middle.
            output[i] = -(betas[i]*S_i*totalInfected) + (sigmas[i-1]*I_im1) - (mu*S_i);
        }

        //Calculate dxdt for each infectious compartment.
        for (unsigned int i=0; i<numStrains; i++) //Note infected compartments are in the second half of the state representation vector.
        {
            const double S_i = currentValues[i];
            const double I_i = currentValues[numStrains+i];

            output[numStrains+i] = +betas[i]*S_i*totalInfected - sigmas[i]*I_i - mu*I_i;
        }

        //Calculate recovered compartment.
        output[numStrains*2] = +sigmas[numStrains-1]*currentValues[(numStrains*2)-1]-mu*R;
        */
    }

    std::vector<double> generate_init_vals(const double _initialInfected)
    {
        std::vector<double> init;

        //Susceptibles
        init.push_back(1.0 - _initialInfected); //Initial naive.
        for (unsigned int i=1; i<numStrains; i++)
            init.push_back(0.0);

        //Infectious
        init.push_back(_initialInfected); //Initial first infections.
        for (unsigned int i=1; i<numStrains; i++)
            init.push_back(0.0);

        //Recovered
        init.push_back(0.0);

        return init;
    }

    void export_num_strains(std::string _name) const
    {
        std::vector<double> temp;
        temp.push_back(numStrains);
        vectorToFile(temp, _name+"_numStrains.csv");
        //std::cout << "Exported numStrains = " << numStrains << "\n";
    }
};
