#pragma once
#include <vector>
#include <string>

//Sweep parameter space - numStrains.
void num_strains_sweep(const unsigned int minStrains, const unsigned int maxStrains, const unsigned int incrementStrains, const double beta0, const double sigma0, const double alpha, const double alphaMax, const double mu, const double initialInfected, const std::string name="si_anon_strains_nocross_numstrains_sweep");
void num_strains_sweep(const std::vector<unsigned int> numStrains, const double beta0, const double sigma0, const double alpha, const double alphaMax, const double mu, const double initialInfected, const std::string name="si_anon_strains_nocross_numstrains_sweep");
