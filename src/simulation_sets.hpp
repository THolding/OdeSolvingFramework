#pragma once
#include <vector>
#include <string>


//Sweep parameter space - numStrains.
std::vector<std::vector<double>> num_strains_sweep(const std::string name, const std::vector<double> beta0s, const double sigma, const double mu, const std::vector<unsigned int> numStrainsList, const double initialInfected, const bool exportAll = false);

std::vector<double> eir_sweep(std::string name, const std::vector<double> beta0s, double sigma, double mu, unsigned int numStrains, double initialInfected,  std::vector<std::vector<double>> &extraData, bool exportAll=false);

//void beta_decay_sweep(const double minBetaDecay, const double maxBetaDecay, const double incrementDecay, const double beta0, const double sigma0, const double mu, const double numCompartments, const double initialInfected);
//void beta_decay_sweep(const std::vector<double> betaDecays, const double beta0, const double sigma0, const double mu, const double numCompartments, const double initialInfected);

void scratchpad();
