#include "simulation_sets.hpp"
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "model_driver.hpp"
#include "si_anon_strains.hpp"
#include "utilities.hpp"


void num_strains_sweep(const unsigned int minStrains, const unsigned int maxStrains, const unsigned int incrementStrains, const double beta0, const double sigma0, const double mu, const double initialInfected)
{
    //Generate numStrains list and call overloaded function.
    std::vector<unsigned int> numStrains;
    for (unsigned int i=minStrains; i<=maxStrains; i+=incrementStrains)
        numStrains.push_back(i);

    num_strains_sweep(numStrains, beta0, sigma0, mu, initialInfected);
}

void num_strains_sweep(const std::vector<unsigned int> numStrains, const double beta0, const double sigma0, const double mu, const double initialInfected)
{
    unsigned int count=0;
    vectorToFile(numStrains, "si_anon_strains_numstrains_list.csv");
    //std::vector<double> prevalence;

    for (unsigned int curNumStrains : numStrains)
    {
        std::cout << "Continuing sweep with numStrains = " << curNumStrains << "...\n";
        SIAnonStrains modelDef(curNumStrains);
        std::vector<double> init = modelDef.generate_init_vals(initialInfected);
        std::vector<double> params = modelDef.generate_param_list(beta0, sigma0, mu);

        ModelDriver model(&modelDef, params, init);
        model.set_dt(0.005);
        model.set_max_time(10000);
        model.set_output_frequency(10);
        model.run("si_anon_strains_numstrains_sweep_"+boost::lexical_cast<std::string>(count));

        model.export_output();

        //Calculate parasite prevalence.
        //std::vector<double> finalVals = model.get_current_values();
        //double tempPrevalence = std::accumulate(finalVals.begin(), finalVals.begin()+curNumStrains, 0);
        //std::cout << "Prevalence = " << tempPrevalence << "\n\n";
        //prevalence.push_back(tempPrevalence);

        count++;
    }

    //vectorToFile(prevalence, "si_anon_strains_numstrains_sweep_prevalence.csv");
}
