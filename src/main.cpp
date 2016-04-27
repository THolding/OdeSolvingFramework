#include <iostream>
#include "model_driver.hpp"
#include "sis_model.hpp"
#include "simulation_sets.hpp"

int main()
{
    ///Example execution of a single simulation.
    /*const double beta = 3.0;
    const double sigma = 1.0;
    const double mu = 0.01; //1.0/(70.0);
    SISModel modelDef;

    std::vector<double> init = {0.999, 0.001, 0.0}; //S, I and R initial values
    std::vector<double> params = {beta, sigma, mu};

    ModelDriver model(&modelDef, params, init);
    model.set_dt(0.005);
    model.set_max_time(1000);
    model.set_output_frequency(10);
    model.run("sis_example");

    model.export_output();

    std::cout << "\n\nBeta = " << beta << "\tSigma = " << sigma << "\tMu = " << mu << "\n";
    std::cout << "\nEquilibrium susceptible = " << model.get_current_values()[0];
    std::cout << "\nEquilibrium prevalence = " << model.get_current_values()[1] << "\n";
    */

    ///Run the EIR vs prevalence simulation (uses beta as substitute for prevalence).
    simulate_eir_vs_prevalence();

    return 0;
}
