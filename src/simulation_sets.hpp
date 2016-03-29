#pragma once
#include <vector>
#include <string>

void run_diversity_sweep(const unsigned int numLoci, const std::vector<unsigned int> numAlleles, const double beta, const double gamma, const double sigma, const double mu, const double initInfected, const std::string name="repertoire_multistrain_diversity_sweep");

void run_random_initial_strain_sweep(const unsigned int numLoci, const unsigned int numAlleles, const std::vector<unsigned int> numStrains, const unsigned int numRepeats, const double beta, const double gamma, const double sigma, const double mu, const double initInfected, const std::string name="repertoire_multistrain_random_initial_strain_sweep");
