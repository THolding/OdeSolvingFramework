#include "simulation_sets.hpp"
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "model_driver.hpp"
#include "sir_anon_strains.hpp"
#include "utilities.hpp"
#include "helper_functions.hpp"

//Returns a list of vectors with:
//  [0] = list of beta0s used in eir sweeps (naive transmission rate).
//  [1...] = one prevalence list for each numStrains used.
std::vector<std::vector<double>> num_strains_sweep(const std::string name, const std::vector<double> beta0s, const double sigma, const double mu, const std::vector<unsigned int> numStrainsList, const double initialInfected, const bool exportAll)
{
    std::vector<std::vector<double>> output;
    std::vector<std::vector<double>> recoveredPops;

    //Must convert numStrainsList from std::vector<unsigned int> to std::vector<double>...
    output.push_back(beta0s);
    for (unsigned int numStrains : numStrainsList)
    {
        std::string curName = name+"_"+boost::lexical_cast<std::string>(numStrains);
        std::cout << "Starting numStrains sweep " << curName << ".\n";

        std::vector<double> prevalences = eir_sweep(curName, beta0s, sigma, mu, numStrains, initialInfected, recoveredPops, exportAll);
        output.push_back(prevalences);
    }

    vectorToFile(numStrainsList, name+"_numStrainsList.csv");

    matrixToFile(recoveredPops, name+"_recovered.csv", ", ");

    return output;
}

//Returns a list of prevalences for each starting beta0 used
std::vector<double> eir_sweep(std::string name, const std::vector<double> beta0s, double sigma, double mu, unsigned int numStrains, double initialInfected,  std::vector<std::vector<double>> &extraData, bool exportAll)
{
    std::cout << "Starting EIR sweep with beta0 = ";

    //Temp data collection
    std::vector<double> recoveredPop;

    std::vector<double> prevalences;
    unsigned int beta0Index = 0;
    for (const double beta0 : beta0s)
    {
        std::cout << beta0 << ", ";

        SIRAnonStrains modelDef(numStrains);

        //Define initial values and parameters
        std::vector<double> init = modelDef.generate_init_vals(initialInfected);
        std::vector<double> betas = linear(beta0, numStrains);
        std::vector<double> sigmas(numStrains);
        std::fill(sigmas.begin(), sigmas.end(), sigma);

        std::vector<double> params; //Place to merge all parameters
        params.insert(params.end(), betas.begin(), betas.end());
        params.insert(params.end(), sigmas.begin(), sigmas.end());
        params.push_back(mu);

        ModelDriver model(&modelDef, params, init);
        model.set_dt(0.005);
        model.set_max_time(10000);
        model.set_output_frequency(10);
        model.run(name);

        //temp calc prev.
        std::vector<double> finalValues = model.get_current_values();
        double prevalence = 0.0;
        for (unsigned int i=numStrains; i<numStrains*2; i++)
            prevalence += finalValues[i];
        recoveredPop.push_back(finalValues.back());
        std::cout << "Calculated prev: " << prevalence << "\tt = " << finalValues[0] << "\n";
        prevalences.push_back(prevalence);

        if (exportAll)
        {
            model.export_output("plots_check/"+name+"_single_"+boost::lexical_cast<std::string>(beta0Index)+"_"+boost::lexical_cast<std::string>(beta0));
            modelDef.export_num_strains(name);
        }
        beta0Index++;
    }

    std::cout << "\nPrevalences:\t";
    for (double d : prevalences) std::cout << d << "\t";
    std::cout << "\n\n";

    extraData.push_back(recoveredPop);

    return prevalences;
}


void scratchpad()
{
    std::string name = "numstrains_sweep";
    std::vector<unsigned int> numStrainsList = {2, 3, 4, 5};//, 8, 12, 20, 30, 50, 100};//{3, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 50, 100};//, 4, 5, 7, 9, 12, 15, 20};
    double beta0 = 1.13;
    const double sigma = 0.25;
    const double mu = 1.0/50.0;
    const double initialInfected = 0.0001;

    for (unsigned int numStrains : numStrainsList)
    {
        SIRAnonStrains modelDef(numStrains);

        std::vector<double> init = modelDef.generate_init_vals(initialInfected);
        std::vector<double> betas = linear(beta0, numStrains);
        std::vector<double> sigmas(numStrains);
        std::fill(sigmas.begin(), sigmas.end(), sigma);

        std::vector<double> params; //Place to merge all parameters
        params.insert(params.end(), betas.begin(), betas.end());
        params.insert(params.end(), sigmas.begin(), sigmas.end());
        params.push_back(mu);

        ModelDriver model(&modelDef, params, init);
        model.set_dt(0.005);
        model.set_max_time(10000);
        model.set_output_frequency(10);
        model.run(name+"_"+boost::lexical_cast<std::string>(numStrains)+"_"+boost::lexical_cast<std::string>(beta0));

        model.export_output();
    }
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
        SIRAnonStrains modelDef(numCompartments);
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
