#include <iostream>
#include "model_driver.hpp"
#include "si_anon_strains.hpp"
#include "simulation_sets.hpp"

int main()
{
    ///Example execution of a single simulation.
    /*const unsigned int numStrains = 10;
    const double beta0 = 2.0;
    const double sigma = 0.5;
    const double mu = 1.0/50.0;
    const double alpha0 = 0.0;
    const double alphaMax = 0.0;
    const std::string runName = "si_anon_strains_nocross";
    SIAnonStrains modelDef(numStrains);

    std::vector<double> init = modelDef.generate_init_vals(0.001);
    std::vector<double> params = modelDef.generate_param_list(beta0, sigma, alpha0, alphaMax, mu);

    modelDef.export_params(params, runName);


    ModelDriver model(&modelDef, params, init);
    model.set_dt(0.005);
    model.set_max_time(1000);
    model.set_output_frequency(10);
    model.run(runName);

    model.export_output();
    modelDef.export_params(params, runName);
    modelDef.export_num_strains(runName);*/

    ///Run parameter sweep.
    const double beta0 = 2.0;
    const double sigma = 0.5;
    const double alpha0 = 0.0;
    const double alphaMax = 0.5;
    const double mu = 1.0/50.0;
    const double initialInfected = 0.0001;
    //const std::vector<unsigned int> numStrains = { 3, 5, 8, 10, 13, 15, 20, 30 };

    num_strains_sweep(5, 200, 3, beta0, sigma, alpha0, 1.0, mu, initialInfected, "si_anon_strains_exp_cross4_numstrains_sweep");
    num_strains_sweep(5, 200, 3, beta0, sigma, alpha0, 0.75, mu, initialInfected, "si_anon_strains_exp_cross3_numstrains_sweep");
    num_strains_sweep(5, 200, 3, beta0, sigma, alpha0, 0.5, mu, initialInfected, "si_anon_strains_exp_cross2_numstrains_sweep");
    num_strains_sweep(5, 200, 3, beta0, sigma, alpha0, 0.25, mu, initialInfected, "si_anon_strains_exp_cross1_numstrains_sweep");
    num_strains_sweep(5, 200, 3, beta0, sigma, alpha0, 0.0, mu, initialInfected, "si_anon_strains_exp_cross0_numstrains_sweep");

    return 0;
}
