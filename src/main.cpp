#include <iostream>
#include "model_driver.hpp"
#include "sis_anon_strains.hpp"
#include "simulation_sets.hpp"

int main()
{
    ///Example execution of a single simulation.
    /*const unsigned int numStrains = 3;
    const double beta0 = 2.0;
    const double sigma = 0.5;
    const double mu = 1.0/50.0;
    const double betaDecayRate = 0.1;
    const double sigmaDecayRate = 0.0;
    SISAnonStrains modelDef(numStrains);

    std::vector<double> init = modelDef.generate_init_vals(0.001);
    std::vector<double> params = modelDef.generate_param_list(beta0, sigma, mu, betaDecayRate, sigmaDecayRate);

    ModelDriver model(&modelDef, params, init);
    model.set_dt(0.005);
    model.set_max_time(1000);
    model.set_output_frequency(10);
    model.run("sis_anon_strains");

    model.export_output();
    modelDef.export_num_strains("sis_anon_strains");*/

    ///Run parameter sweep.
    const double beta0 = 2.0;
    const double sigma = 0.5;
    const double mu = 1.0/50.0;
    const double initialInfected = 0.0001;
    const unsigned int numCompartments = 50;
    const std::vector<double> decayRates = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    beta_decay_sweep(0.0, 2.0, 0.025, beta0, sigma, mu, numCompartments, initialInfected);

    return 0;
}
