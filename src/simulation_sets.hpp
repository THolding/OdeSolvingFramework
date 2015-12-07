#pragma once
#include <vector>


//Sweep parameter space - numStrains.
void num_strains_sweep(const unsigned int minStrains, const unsigned int maxStrains, const unsigned int incrementStrains, const double beta0, const double sigma0, const double mu, const double initialInfected);
void num_strains_sweep(const std::vector<unsigned int> numStrains, const double beta0, const double sigma0, const double mu, const double initialInfected);
