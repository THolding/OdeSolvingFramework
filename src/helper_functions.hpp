#pragma once
#include <vector>

std::vector<double> linear(const double beta0, const unsigned int numStrains);

std::vector<double> calc_equal_range(const double startpoint, const double interval, const double endpoint);
