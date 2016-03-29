#pragma once

//Example ross-macdonald model
class ExampleModel : public ModelDefinition
{
public:
    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        float H = params[0]; //host population size
        float M = params[1]; //mosquito population size
        float a = params[2]; //Mosquito bite rate to humans
        float b = params[3]; //probability infectious bite leads to infection in humans
        float c = params[4]; //probability that a mosquito biting infectious human becomes infectious
        float mu_m = params[5]; //Mosquito death rate (1-p)
        float mu_h = params[6]; //host death rate
        float tau = params[7]; //rate infected mosquitoes become infectious
        float sigma = params[8]; //recovery rate of hosts.

        float m = M/H;
        float p = 1.0 - mu_m;
        float n = 1.0/tau;

        double Z = currentValues[0]; //Infectious hosts
        double X = currentValues[1]; //Latent vectors
        double Y = currentValues[2]; //Infectious vectors

        //Infectious hosts... (Z)
        output[0] = a*c*Y*(1.0 - Z) - sigma*Z - mu_h*Z;

        //Infected (latent) vectors (X)
        output[1] = m*a*b*Z*(1.0 - Y - X) - tau*X - mu_m*X;

        //Infectious vectors (Y)
        output[2] = tau*X - mu_m*Y;
    }
};
