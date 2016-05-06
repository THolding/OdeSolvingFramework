#pragma once
#include <vector>
#include <string>


//Using an array of beta0's, will initialise and run the model for each beta0 and output the equilibrium prevalence at each beta 0 to file.
std::vector<double> eir_sweep(std::string name, const std::vector<double> beta0s, double sigma, double mu, unsigned int numStrains, double initialInfected,  std::vector<std::vector<double>> &extraData, bool exportAll=false);

//Sweep parameter space - numStrains.
//Repeats eir_sweep for each numStrains value in numStrainsList.
std::vector<std::vector<double>> eir_sweep_each_numstrains(const std::string name, const std::vector<double> beta0s, const double sigma, const double mu, const std::vector<unsigned int> numStrainsList, const double initialInfected, const bool exportAll = false);


//Using an array of numStrains, will initialise and run the model for each numStrains and output the equilibrium prevalence at each beta 0 to file.
std::vector<double> num_strains_sweep(std::string name, const std::vector<unsigned int> numStrainsList, double beta0, double sigma, double mu, double initialInfected, bool exportAll, bool suppressOutput);

//Sweep parameter space - beta0
//Repeats num_strains_sweep for each beta0 value in beta0List.
std::vector<std::vector<double>> num_strains_sweep_each_beta0(const std::string name, const std::vector<unsigned int> numStrainsList, const std::vector<double> beta0List, const double sigma, const double mu, const double initialInfected, const bool exportAll = false, const bool suppressOutput = false);


void scratchpad();
