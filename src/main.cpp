#include <iostream>
#include "model_driver.hpp"
#include "sir_model.hpp"
#include "simulation_sets.hpp"

int main()
{
    ///Example execution of a single simulation.
    /*
    const double beta = 3.0;
    const double sigma = 1.0;
    const double mu = 0.01; //1.0/(70.0);
    SIRModel modelDef;

    std::vector<double> init = {0.999, 0.001, 0.0}; //S, I and R initial values
    std::vector<double> params = {beta, sigma, mu};

    ModelDriver model(&modelDef, params, init);
    model.set_dt(0.005);
    model.set_max_time(1000);
    model.set_output_frequency(10);
    model.run("sir_example");

    model.export_output();

    std::cout << "\n\nBeta = " << beta << "\tSigma = " << sigma << "\tMu = " << mu << "\n";
    std::cout << "\nEquilibrium susceptible = " << model.get_current_values()[0];
    std::cout << "\nEquilibrium prevalence = " << model.get_current_values()[1];
    std::cout << "\nEquilibrium recovered = " << model.get_current_values()[2] << "\n";
    */


    ///Run the EIR vs prevalence simulation (uses beta as substitute for prevalence).
    simulate_eir_vs_prevalence();

    ///Run parameter sweep.
    /*const double beta0 = 2.0;
    const double sigma = 0.5;
    const double mu = 1.0/50.0;
    const double initialInfected = 0.0001;
    const unsigned int numCompartments = 50;
    const std::vector<double> decayRates = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    beta_decay_sweep(0.0, 2.0, 0.025, beta0, sigma, mu, numCompartments, initialInfected);
    */

    return 0;
}
