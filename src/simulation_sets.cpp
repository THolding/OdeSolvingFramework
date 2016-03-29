#include "simulation_sets.hpp"
#include <random>
#include <ctime>
#include "repertoire_multistrain_w.hpp"
#include "model_driver.hpp"
#include <boost/lexical_cast.hpp>

void run_diversity_sweep(const unsigned int numLoci, const std::vector<unsigned int> numAlleles, const double beta, const double gamma, const double sigma, const double mu, const double initInfected, const std::string name)
{
    std::srand(std::time(NULL));

    unsigned int count=0;
    for (const unsigned int curNumAlleles : numAlleles) //For each level of diversity we're interested in...
    {
        RepertoireMultistrainW modelDef(numLoci, curNumAlleles);

        //Initialise parameters/
        std::vector<double> params;
        modelDef.calculate_parameter_set(beta, gamma, sigma, mu, params);
        std::cout << "params:\n";
        for (unsigned int i=0; i<params.size(); i++)
            std::cout << i << ":\t" << params[i] << "\n";
        std::cout << "\n\n\n";

        //Create the initial values for.
        const double initPerStrain = initInfected/(double)modelDef.get_num_strains();
        std::vector<double> ys;
        for (unsigned int i=0; i<modelDef.get_num_strains(); i++)
            ys.push_back(initPerStrain);
        for (unsigned int i=0; i<ys.size(); i++)
        {
            if (randDouble() < 0.5)
            {
                ys[i] += (initPerStrain*0.001); //Add a little
                ys[rand()%modelDef.get_num_strains()] -= (initPerStrain*0.001);//Remove a little from another strain to keep total the same.
            }
        }

        //How use initial infectious sizes and params to calculate the full set of initial conditions.
        std::vector<double> init;
        modelDef.calculate_initial_values(ys, params, init);

        //Debugging: output initial values.
        std::cout << "Initial values:\n";
        for (unsigned int i=0; i<init.size(); i++)
            std::cout << i << ":\t" << init[i] << "\n";

        ModelDriver model(&modelDef, params, init);
        model.set_dt(0.05);
        model.set_max_time(50000);
        model.set_output_frequency(20);
        //model.set_stop_threshold(0);
        model.run(name + boost::lexical_cast<std::string>(count));
        modelDef.export_num_strains(name + boost::lexical_cast<std::string>(count));

        model.export_output();
        count++;
    }
}

void run_random_initial_strain_sweep(const unsigned int numLoci, const unsigned int numAlleles, const std::vector<unsigned int> numStrains, const unsigned int numRepeats, const double beta, const double gamma, const double sigma, const double mu, const double initInfected, const std::string name)
{
    std::srand(std::time(NULL));

    std::vector<double> prevalences;

    unsigned int count=0;
    for (const unsigned int curNumStrains : numStrains) //For each level of diversity we're interested in...
    {
        double prevalence = 0.0; //Prevalence for this number of strains (to be averaged over all reps.
        for (unsigned int rep=0; rep<numRepeats; rep++)
        {
            RepertoireMultistrainW modelDef(numLoci, numAlleles);

            if (curNumStrains > modelDef.get_num_strains())
                throw std::runtime_error("Cannot initialise more strains than possible antigen combinations in run_random_initial_strain_sweep().\n");

            //Initialise parameters
            std::vector<double> params;
            modelDef.calculate_parameter_set(beta, gamma, sigma, mu, params);
            //std::cout << "params:\n";
            //for (unsigned int i=0; i<params.size(); i++)
            //    std::cout << i << ":\t" << params[i] << "\n";
           // std::cout << "\n\n\n";

            //Create the initial values for.
            const double initPerStrain = initInfected/(double)curNumStrains; //Initial infectious for each non-zero strain.

            //Get list of strains to select for initialisation.
            std::vector<unsigned int> strainsToInitialise;
            for (unsigned int i=0; i<curNumStrains; i++)
                strainsToInitialise.push_back(randomNumberWithExclusions(modelDef.get_num_strains(), strainsToInitialise));

            //Initialise ys with all zeros.
            std::vector<double> ys;
            for (unsigned int i=0; i<modelDef.get_num_strains(); i++)
                ys.push_back(0.0);

            //Initialise selected random strains.
            for (const unsigned int i : strainsToInitialise)
            {
                double noise = 0.001 * ((randDouble()*2)-1);
                ys[i] = initPerStrain + (initPerStrain*(noise)); //Add a 0.1% random noise.
            }

            //How use initial infectious sizes and params to calculate the full set of initial conditions.
            std::vector<double> init;
            modelDef.calculate_initial_values(ys, params, init);

            //Debugging: output initial values.
            //std::cout << "Initial values:\n";
            //for (unsigned int i=0; i<init.size(); i++)
            //    std::cout << i << ":\t" << init[i] << "\n";

            ModelDriver model(&modelDef, params, init);
            model.set_dt(0.1);
            model.set_max_time(10000);
            model.set_output_frequency(30);
            //model.set_stop_threshold(0);
            model.run(name + boost::lexical_cast<std::string>(count));
            //modelDef.export_num_strains(name + boost::lexical_cast<std::string>(count));

            std::vector<double> finalValues = model.get_current_values();
            double tempPrevalence = 0.0;
            for (unsigned int i=modelDef.get_num_strains()*2; i<modelDef.get_num_strains()*3; i++) //Sum over all infectious compartments.
                tempPrevalence += finalValues[i];
            prevalence += tempPrevalence;

            //model.export_output();
            std::cout << "End num strains: " << curNumStrains << " rep: " << rep << "\n";
        } //End rep

        //Calculate mean prevalence.
        prevalence = prevalence / (double) numRepeats;
        prevalences.push_back(prevalence);

        count++;
    }

    vectorToFile(numStrains, name+"_selectedStrains.csv");
    vectorToFile(prevalences, name+"_prevalence.csv");
}
