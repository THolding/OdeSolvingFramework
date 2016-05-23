#pragma once
#include <iostream>
#include <cmath>
#include "utilities.hpp"
#include "model_definition.hpp"

//This is a standard two strain model system with which to compare the SIRAnonStrains with to ensure equivalency.
class MultistrainEquivalent : public ModelDefinition
{
private:
public:
    MultistrainEquivalent() {  }
    //Parameters: beta1, beta2, sigma1, sigma2, mu.
    //Compartments as follows: S, S1, S2, I1, I2, I12, I21, R.
    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        //Define parameter constants for easier readability.
        const double S   = currentValues[0]; //Naive individuals
        const double S1  = currentValues[1];
        const double S2  = currentValues[2];
        const double I1  = currentValues[3];
        const double I2  = currentValues[4];
        const double I12 = currentValues[5];
        const double I21 = currentValues[6];
        const double R   = currentValues[7];

        const double beta1  = params[0];
        const double beta2  = params[1];
        const double sigma1 = params[2];
        const double sigma2 = params[3];
        const double mu     = params[4];


        //Calculate dxdt for naive individuals (S).
        output[0] = mu - beta1*S*(I1+I21) - mu*S; //S
        output[1] = sigma1*I1 - beta2*S1*(I2+I12) - mu*S1; //S1
        output[2] = sigma2*I2 - beta1*S2*(I1+I21) - mu*S2; //S2

        output[3] = beta1*S*(I1+I21) - sigma1*I1 - mu*I1; //I1
        output[4] = beta2*S*(I2+I12) - sigma2*I2 - mu*I2; //I2
        output[5] = beta2*S2*(I2+I12) - sigma2*I12 - mu*I12; //I12
        output[6] = beta1*S2*(I1+I21) - sigma1*I21 - mu*I21; //I21

        output[7] = sigma1*I21 + sigma2*I12 - mu*R; //R
    }

    std::vector<double> generate_init_vals(const double _initialInfected)
    {
        double eachStrainInit = _initialInfected/2.0;

        std::vector<double> init;
        init.push_back(1.0-(2*eachStrainInit)); //S
        init.push_back(0.0); //S1
        init.push_back(0.0); //S2
        init.push_back(eachStrainInit); //I1
        init.push_back(eachStrainInit); //I2
        init.push_back(0.0); //I12
        init.push_back(0.0); //I21
        init.push_back(0.0); //IR

        return init;
    }
};

void run_example_multistrain_equivalent()
{
    const double beta1 = 1.0;
    const double beta2 = 1.0;
    const double sigma1 = 0.25;
    const double sigma2 = 0.25;
    const double mu = 1.0/50.0;
    const double initialInfected = 0.001;

    MultistrainEquivalent modelDef;

    //Define initial values and parameters
    std::vector<double> init = modelDef.generate_init_vals(initialInfected); // {0.999, 0.001, 0, 0, 0};
    std::cout << "initial values: ";
    for (double d : init) std::cout << d << " ";
    std::cout << "\n";

    std::vector<double> params = {beta1, beta2, sigma1, sigma2, mu};

    std::cout << "\n\nparams:\n";
    for (double d : params) std::cout << d << " ";
    std::cout << "\n\n\n";

    ModelDriver model(&modelDef, params, init);
    model.set_dt(0.005);
    model.set_max_time(1000);
    model.set_output_frequency(10);
    model.run("multistrain_equivalence");

    //model.export_output();

    //calc prev and output final values
    std::vector<double> finalValues = model.get_current_values();

    //double prevalence = 0;
    //for (unsigned int i=numStrains; i<numStrains*2; i++)
    //    prevalence += finalValues[i];
    //std::cout << "calculated prev: " << prevalence << "\n\n";

    std::cout << "Final Values: ";
    for (double d : finalValues) std::cout << d << " ";
    std::cout << "\n\n";
}
