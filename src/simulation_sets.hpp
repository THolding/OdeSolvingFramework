#pragma once
#include <vector>


//Generates data for an EIR vs prevalence plot, with beta acting as a proxy for EIR.
void simulate_eir_vs_prevalence();


///Obscelete
//Sweep parameter space - numStrains.
//void beta_decay_sweep(const double minBetaDecay, const double maxBetaDecay, const double incrementDecay, const double beta0, const double sigma0, const double mu, const double numCompartments, const double initialInfected);
//void beta_decay_sweep(const std::vector<double> betaDecays, const double beta0, const double sigma0, const double mu, const double numCompartments, const double initialInfected);
