#include "helper_functions.hpp"
#include <vector>

std::vector<double> linear(const double beta0, const unsigned int numStrains)
{
    std::vector<double> betas;
    for (unsigned int i=0; i<numStrains; i++)
        betas.push_back( ((double)(numStrains-i))/((double)numStrains) * beta0 );

    return betas;
}

std::vector<double> calc_equal_range(const double startpoint, const double interval, const double endpoint)
{
    std::vector<double> output;
    for (double d = startpoint; d<=endpoint; d+=interval)
    {
        output.push_back(d);
    }
    return output;
}
