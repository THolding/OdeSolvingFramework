#pragma once
#include <vector>

std::vector<double> linear(const double beta0, const unsigned int numStrains); //Variant specific immunity...

std::vector<double> crossimmunity_exponential(const double beta0, const unsigned int numStrains, const double crossreactivity); //Exponential transmission to model cross-reactivity plus variant specific immunity.

std::vector<double> calc_equal_range(const double startpoint, const double interval, const double endpoint);

std::vector<double> calc_recovery_increase(const double naiveSigma, const unsigned int numStrains, const double maxSigma, const double rate);
