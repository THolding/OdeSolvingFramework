#pragma once
#include <algorithm>
#include "utilities.hpp"

class RepertoireMultistrainW : public ModelDefinition
{
private:
    const unsigned int repertoireSize;
    const unsigned int numVariants;
    unsigned int numStrains;
    std::vector<std::vector<unsigned int>> overlapMatrix;
    std::vector<unsigned int> antigens;
    std::vector<std::vector<unsigned int>> strains;
public:
    RepertoireMultistrainW(unsigned int _repertoireSize, unsigned int _numVariants) : repertoireSize(_repertoireSize), numVariants(_numVariants)
    {
        generate_antigens();
        generate_strains();
        calculate_overlap_matrix();
    }

    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        const std::vector<double> betas(params.begin(), params.end()-3);
        const double gamma = params.at(params.size()-3);
        const double sigma = params.at(params.size()-2);
        const double mu = params.at(params.size()-1);

        std::cout << "CurVal: ";
        for (auto i:currentValues)
            std::cout << i << ", ";
        std::cout << "\n\n";

        for (std::size_t i=0; i<numStrains; i++)
        {
            const double z_i = currentValues[i];
            const double w_i = currentValues[(numStrains)+i];
            const double y_i = currentValues[(numStrains*2)+i];
            const double foi = betas[i]*y_i; //Force of infection = beta*y_i;

            //dz_i/dt = hosts immune to strain i.
            output[i] = (1.0 - z_i)*foi - mu*z_i;

            //dw_i/dt = hosts immune to strains overlapping with i (including i).
            output[numStrains+i] = (1.0 - w_i)*phi(i, currentValues, betas) - mu*w_i;

            //dy_i/dt = hosts infectious with strain i.
            output[(numStrains*2)+i] = ((1.0 - w_i) + (1.0 - gamma)*(w_i - z_i))*foi - sigma*y_i;
        }
    }

    void calculate_overlap_matrix()
    {
        //Calculate overlap matrix.
        for (std::size_t curStrain=0; curStrain < numStrains; curStrain++) //For each strain store a list of overlap sizes corresponding to each other strain.
        {
            std::vector<unsigned int> row;
            for (std::size_t overlapStrain=0; overlapStrain < strains.size(); overlapStrain++) //For each strain, calculate the degree of overlap to the current strain.
            {
                //No cross-reactive immunity to self, since self is subject to direct immunity.
                //if (curStrain == overlapStrain)
                //{
                //    row.push_back(0);
                //    continue;
                //}

                unsigned int overlapSize = 0;
                for (unsigned int antigen : strains[overlapStrain]) //Search for each antigen from the overlapStrain in the current strain and count overlaps.
                {
                    if (std::find(strains[curStrain].begin(), strains[curStrain].end(), antigen) != strains[curStrain].end())
                        overlapSize++;
                }
                row.push_back(overlapSize);
            }
            overlapMatrix.push_back(row);
        }
    }

    double phi(const unsigned int _strain, const std::vector<double>& _betas, const std::vector<double>& _curVals) const
    {
        double output = 0.0;
        for (std::size_t j=0; j<numStrains; j++)
        {
            if (overlapMatrix[_strain][j] > 0)
                output += _betas[j]*_curVals[(numStrains*2)+j]; //+= force of infection.
        }

        return output;
    }

    void generate_antigens()
    {
        for (unsigned int i=0; i<numVariants; i++)
            antigens.push_back(i);
    }

    void generate_strains()
    {
        //Generate list of all possible strains.
        recursive_combinations(repertoireSize, antigens, strains);

        numStrains = strains.size();
    }

    void recursive_combinations(unsigned int _repertoireSize, const std::vector<unsigned int>& _antigens, std::vector<std::vector<unsigned int>>& _strains, std::vector<unsigned int> _curCombination = std::vector<unsigned int>(), unsigned int _offset=0)//, const std::vector<unsigned int> _antigens, std::vector<std::vector<unsigned int>> _strains)
    {
        if (_repertoireSize == 0)
        {
            _strains.push_back(_curCombination);
            return;
        }
        for (unsigned int i = _offset; i <= _antigens.size() - _repertoireSize; ++i)
        {
            _curCombination.push_back(_antigens.at(i));
            recursive_combinations(_repertoireSize-1, _antigens, _strains, _curCombination, i+1);
            _curCombination.pop_back();
        }
    }

    unsigned int get_num_strains() const
    {
        return numStrains;
    }

    void export_num_strains(std::string _name) const
    {
        std::vector<double> temp;
        temp.push_back(numStrains);
        vectorToFile(temp, _name+"_numStrains.csv");
        std::cout << "Exported numStrains = " << get_num_strains() << "\n";
    }
};
