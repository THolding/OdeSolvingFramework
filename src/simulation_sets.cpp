#include "simulation_sets.hpp"
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "model_driver.hpp"
#include "sir_model.hpp"
#include "utilities.hpp"



//Generates data for an EIR vs prevalence plot, with beta acting as a proxy for EIR.
void simulate_eir_vs_prevalence()
{
    //Parameter initialisation
    std::vector<double> init = {0.999, 0.001, 0.0};

    const double sigma = 0.05; //1.0/120.0;
    const double mu = 0.0025; //1.0/90.0;
    std::vector<double> betas;
    for (double i=0.05; i<=2.5; i+=0.05)
        betas.push_back(i);

    //Run the model for each beta value and extract equilibrium prevalence.
    std::vector<double> prevalences;
    for (unsigned int index=0; index<betas.size(); index++)
    {
        //Create and run model for current beta
        SIRModel modelDef;
        std::vector<double> params = {betas[index], sigma, mu};
        ModelDriver model(&modelDef, params, init);
        model.set_dt(0.005);
        model.set_max_time(10000);
        model.set_output_frequency(1);
        model.run("sir_eirvsprev_"+boost::lexical_cast<std::string>(index));
        model.export_output();

        //Calculate prevalence
        std::vector<double> finalVals = model.get_current_values();
        prevalences.push_back(finalVals[2]);
    }

    //Output betas vs prevalence
    std::vector<std::vector<double>> output;
    output.push_back(betas);
    output.push_back(prevalences);
    matrixToFile(output, "eirvsprev_prevalences.csv", ", ");
}








/*void beta_decay_sweep(const double minBetaDecay, const double maxBetaDecay, const double incrementDecay, const double beta0, const double sigma0, const double mu, const double numCompartments, const double initialInfected)
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
}*/
