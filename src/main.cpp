#include <iostream>
#include "model_driver.hpp"
#include "sir_anon_strains.hpp"
#include "simulation_sets.hpp"
#include "helper_functions.hpp"

int main()
{

    ///Example execution of a single simulation.
    /*const unsigned int numStrains = 2;
    const double beta0 = 1.0;
    const double sigma = 0.25;
    const double mu = 1.0/50.0;
    const double initialInfected = 0.001;
    //const double betaDecayRate = 0.1;
    //const double sigmaDecayRate = 0.0;
    SIRAnonStrains modelDef(numStrains);

    //Define initial values and parameters
    std::vector<double> init = modelDef.generate_init_vals(initialInfected); // {0.999, 0.001, 0, 0, 0};
    std::cout << "initial values: ";
    for (double d : init) std::cout << d << " ";
    std::cout << "\n";

    std::vector<double> betas = linear(beta0, numStrains);
    std::vector<double> sigmas(numStrains);
    std::fill(sigmas.begin(), sigmas.end(), sigma);

    std::vector<double> params; //Place to merge all parameters
    params.insert(params.end(), betas.begin(), betas.end());
    params.insert(params.end(), sigmas.begin(), sigmas.end());
    params.push_back(mu);

    std::cout << "\n\nparams:\n";
    for (double d : params) std::cout << d << " ";
    std::cout << "\n";

    ModelDriver model(&modelDef, params, init);
    model.set_dt(0.005);
    model.set_max_time(1000);
    model.set_output_frequency(10);
    model.run("sir_anon_strains");

    model.export_output();
    modelDef.export_num_strains("sir_anon_strains");

    //temp calc prev.
    std::vector<double> finalValues = model.get_current_values();
    //std::cout << "finalValues:\n";
    //for (double d : finalValues) std::cout << d << "\t";
    //std::cout << "\n\n";
    double prevalence = 0;
    for (unsigned int i=numStrains; i<numStrains*2; i++)
        prevalence += finalValues[i];
    std::cout << "calculated prev: " << prevalence << "\n\n";
    std::cout << "Final Values: ";
    for (double d : finalValues) std::cout << d << " ";
    std::cout << "\n\n";
    */

    std::vector<double> temp = calc_recovery_increase(1, 10, 5, 0.5);
    for (double d : temp) std::cout << d << "\n";

    ///Run a series of eir sweeps with different numstrains
    //scratchpad();

    std::string name = "numstrains_sweep_crossimmunity";
    std::vector<unsigned int> numStrainsList = {2, 4, 6, 8, 12, 20, 30, 50, 100};//{3, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 50, 100};//, 4, 5, 7, 9, 12, 15, 20};
    std::vector<double> beta0s = calc_equal_range(0.2, 0.05, 3.0); //Naive transmission rates to sweep through.
    const double sigma = 0.25;
    const double mu = 1.0/50.0;
    const double initialInfected = 0.0001;
    std::vector<std::vector<double>> output = num_strains_sweep(name, beta0s, sigma, mu, numStrainsList, initialInfected, false);

    matrixToFile(output, name+".csv", ", ");

    return 0;
}
