#include "simulation_sets.hpp"
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "model_driver.hpp"
#include "si_anon_strains.hpp"
#include "utilities.hpp"


void num_strains_sweep(const unsigned int minStrains, const unsigned int maxStrains, const unsigned int incrementStrains, const double beta0, const double sigma0, const double alpha, const double alphaMax, const double mu, const double initialInfected, const std::string name)
{
    //Generate numStrains list and call overloaded function.
    std::vector<unsigned int> numStrains;
    for (unsigned int i=minStrains; i<=maxStrains; i+=incrementStrains)
        numStrains.push_back(i);

    num_strains_sweep(numStrains, beta0, sigma0, alpha, alphaMax, mu, initialInfected, name);
}

void num_strains_sweep(const std::vector<unsigned int> numStrains, const double beta0, const double sigma0, const double alpha, const double alphaMax, const double mu, const double initialInfected, const std::string name)
{
    unsigned int count=0;
    vectorToFile(numStrains, "si_anon_strains_numstrains_list.csv");
    std::vector<double> prevalence;

    for (unsigned int curNumStrains : numStrains)
    {
        std::cout << "Continuing sweep with " << name << " = " << curNumStrains << "...\n";
        SIAnonStrains modelDef(curNumStrains);
        std::vector<double> init = modelDef.generate_init_vals(initialInfected);
        std::vector<double> params = modelDef.generate_param_list(beta0, sigma0, alpha, alphaMax, mu);

        ModelDriver model(&modelDef, params, init);
        model.set_dt(0.005);
        model.set_max_time(10000);
        model.set_output_frequency(10);
        model.run(name+"_"+boost::lexical_cast<std::string>(count));

        if (false)
            model.export_output();

        if (curNumStrains == numStrains[0]) //If first iteration...
            modelDef.export_params(params, name);

        //Calculate parasite prevalence.
        std::vector<double> finalVals = model.get_current_values();
        double tempPrevalence = 0.0;
        for (unsigned int i=curNumStrains+1; i<(curNumStrains*2)+1; i++)
            tempPrevalence+=finalVals[i];
        std::cout << "\nPrevalence = " << tempPrevalence << "\n";
        prevalence.push_back(tempPrevalence);

        //Calculate recovered...
        //recovered.push_back(finalVals[curNumStrains]);

        count++;
    }

    vectorToFile(prevalence, name+"_prevalence.csv");
    //vectorToFile(recovered, name+"_recovered.csv");
}
