#include "simulation_sets.hpp"
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "model_driver.hpp"
#include "sis_anon_strains.hpp"
#include "utilities.hpp"


void beta_decay_sweep(const double minBetaDecay, const double maxBetaDecay, const double incrementDecay, const double beta0, const double sigma0, const double mu, const double numCompartments, const double initialInfected)
{
    //Generate numStrains list and call overloaded function.
    std::vector<double> betaDecays;
    for (double i=minBetaDecay; i<=maxBetaDecay; i+=incrementDecay)
        betaDecays.push_back(i);

    beta_decay_sweep(betaDecays, beta0, sigma0, mu, numCompartments, initialInfected);
}

void beta_decay_sweep(const std::vector<double> betaDecays, const double beta0, const double sigma0, const double mu, const double numCompartments, const double initialInfected)
{
    unsigned int count=0;
    vectorToFile(betaDecays, "sis_anon_strains_betadecays_list.csv");
    //std::vector<double> prevalence;

    for (double curDecay : betaDecays)
    {
        std::cout << "Continuing sweep with betaDecay = " << curDecay << "...\n";
        SISAnonStrains modelDef(numCompartments);
        std::vector<double> init = modelDef.generate_init_vals(initialInfected);
        std::vector<double> params = modelDef.generate_param_list(beta0, sigma0, mu, curDecay, 0.0);

        ModelDriver model(&modelDef, params, init);
        model.set_dt(0.005);
        model.set_max_time(10000);
        model.set_output_frequency(10);
        model.run("sis_anon_strains_betadecay_sweep_"+boost::lexical_cast<std::string>(count));

        model.export_output();
        modelDef.export_num_strains("sis_anon_strains");

        //Calculate parasite prevalence.
        //std::vector<double> finalVals = model.get_current_values();
        //double tempPrevalence = std::accumulate(finalVals.begin(), finalVals.begin()+curDecay, 0);
        //std::cout << "Prevalence = " << tempPrevalence << "\n\n";
        //prevalence.push_back(tempPrevalence);

        count++;
    }

    //vectorToFile(prevalence, "si_anon_strains_numstrains_sweep_prevalence.csv");
}
