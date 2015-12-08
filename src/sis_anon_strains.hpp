#pragma once
#include <iostream>
#include <cmath>
#include "utilities.hpp"
#include "model_definition.hpp"

class SISAnonStrains : public ModelDefinition
{
private:
    const unsigned int numStrains;
public:
    SISAnonStrains(const unsigned int _numStrains) : numStrains(_numStrains)
    {  }


    //Parameters: numStrains+1 betas (final beta is typically 0, indicating total immunity).
    //            numStrains sigmas.
    //            mu at params.end()-1.
    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        const std::vector<double> betas(params.begin(), params.begin()+numStrains+1);
        const std::vector<double> sigmas(params.begin()+numStrains+1, params.begin()+(2*numStrains)+1);
        const double mu = params.back();
        const double N = 1.0;

        double totalInfected = 0;
        for (unsigned int i=0; i<numStrains; i++)
            totalInfected += currentValues[i+numStrains];

        ////Calculate dxdt for each compartment.
        //Calculate dxdt for the first susceptible compartment.
        output[0] = -betas[0]*currentValues[0]*totalInfected +mu*N -mu*currentValues[0]; //First compartment is a special case because it has no inflow flux from an infected compartment.

        //Calculate dxdt for middle susceptible compartments.
        for (unsigned int i=1; i<numStrains-1; i++)
        {
            const double S_i = currentValues[i];
            const double I_im1 = currentValues[numStrains+i-1];

            //Note that the first 50% of compartments are taken to be the susceptible compartments, while the second 50% are infected.
            output[i] = -(betas[i]*S_i*totalInfected) + (sigmas[i-1]*I_im1) - (mu*S_i);
        }

        //Calculate final susceptible compartment.
        output[numStrains-1] = -(betas[numStrains-1]*currentValues[numStrains-1]*totalInfected) + (sigmas[numStrains-2]*currentValues[(numStrains*2)-2]) + (sigmas[numStrains-1]*currentValues[(numStrains*2)-1]) - (mu*currentValues[numStrains-1]);

        //Calculate dxdt for each infectious compartment.
        for (unsigned int i=0; i<numStrains; i++) //Note infected compartments are in the second half of the state representation vector.
        {
            const double S_i = currentValues[i];
            const double I_i = currentValues[numStrains+i];

            output[numStrains+i] = +betas[i]*S_i*totalInfected - sigmas[i]*I_i - mu*I_i;
        }
    }

    std::vector<double> generate_param_list(const double _beta0, const double _sigma0, const double _mu, const double _betaDecay, const double _sigmaDecay = 0.0)
    {
        std::vector<double> params;
        //Add beta_i parameters
        for (unsigned int i=0; i<=numStrains; i++)
            params.push_back(_beta0 * std::exp(-_betaDecay*(double)i));

        //Add sigma parameters
        for (unsigned int i=0; i<=numStrains; i++)
            params.push_back(_sigma0 * std::exp(-_sigmaDecay*(double)i));

        //Add mu parameter
        params.push_back(_mu);

        return params;
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

        return init;
    }

    void export_num_strains(std::string _name) const
    {
        std::vector<double> temp;
        temp.push_back(numStrains);
        vectorToFile(temp, _name+"_numStrains.csv");
        std::cout << "Exported numStrains = " << numStrains << "\n";
    }
};
