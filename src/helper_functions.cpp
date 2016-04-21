#include "helper_functions.hpp"
#include <vector>
#include <cmath>

std::vector<double> linear(const double beta0, const unsigned int numStrains)
{
    std::vector<double> betas;
    for (unsigned int i=0; i<numStrains; i++)
        betas.push_back( ((double)(numStrains-i))/((double)numStrains) * beta0 );

    return betas;
}

//Exponential cross-protection plus variant specific immunity... (transmission rate).
std::vector<double> crossimmunity_exponential(const double beta0, const unsigned int numStrains, const double crossreactivity)
{
    std::vector<double> betas;
    for (unsigned int i=0; i<numStrains; i++)
        betas.push_back( ((double)(numStrains-i))/((double)numStrains) * beta0 * std::exp(-(double)i * crossreactivity));

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

std::vector<double> calc_recovery_increase(const double naiveSigma, const unsigned int numStrains, const double maxSigma, const double rate)
{
    std::vector<double> sigmas;
    for (unsigned int i=0; i<numStrains; i++)
        sigmas.push_back(naiveSigma + (1.0 - std::exp(-(double)i * rate))*maxSigma );

    return sigmas;
}
