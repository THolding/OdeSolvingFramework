#pragma once
#include <cmath>

class GuptaMultistrainW : public ModelDefinition
{
private:
    unsigned int numLoci;
    unsigned int numAlleles;
    std::vector<std::vector<unsigned int>> overlapMatrix;
    unsigned int numStrains;
public:
    GuptaMultistrainW(unsigned int _numLoci, unsigned int _numAlleles) :
        numLoci(_numLoci), numAlleles(_numAlleles), numStrains(std::pow(_numLoci, _numAlleles))
    {
        calculate_overlap_matrix();
    }

    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        const std::vector<double> betas(params.begin(), params.end()-3);
        const double gamma = params.at(params.size()-3);
        const double sigma = params.at(params.size()-2);
        const double mu = params.at(params.size()-1);

        /*
        //Lamba indices.
        const unsigned int AX = 0;
        const unsigned int AY = 1;
        const unsigned int BX = 2;
        const unsigned int BY = 3;

        //Proportion immune to strain i:
        const double z_ax = currentValues[0]; //Proportion hosts immune to ax.
        const double z_ay = currentValues[1]; //Proportion hosts immune to ay.
        const double z_bx = currentValues[2]; //Proportion hosts immune to bx.
        const double z_by = currentValues[3]; //Proportion hosts immune to by.

        //Proportions of hosts immune to any strain sharing alleles with i (including i).
        const double w_ax = currentValues[4]; //Proportion hosts immune to strains overlapping ax (including self).
        const double w_ay = currentValues[5]; //Proportion hosts immune to strains overlapping ay (including self).
        const double w_bx = currentValues[6]; //Proportion hosts immune to strains overlapping bx (including self).
        const double w_by = currentValues[7]; //Proportion hosts immune to strains overlapping by (including self).

        //Proportions of hosts infectious with strain i.
        const double y_ax = currentValues[8]; //Proportion hosts infectious with strain ax.
        const double y_ay = currentValues[9]; //Proportion hosts infectious with strain ay.
        const double y_bx = currentValues[10]; //Proportion hosts infectious with strain bx.
        const double y_by = currentValues[11]; //Proportion hosts infectious with strain by.

        //Calculate force of infections.
        const std::vector<double> lamba = { beta[0]*y_ax, beta[1]*y_ay, beta[2]*y_bx, beta[3]*y_by };
        */
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

        /*
        //dz/dt: Change to proportion of hosts immune to strains.
        output[0] = (1.0 - z_ax)*lamba[AX] - mu*z_ax; //z_ax/dt.
        output[1] = (1.0 - z_ay)*lamba[AY] - mu*z_ay; //z_ay/dt.
        output[2] = (1.0 - z_bx)*lamba[BX] - mu*z_bx; //z_bx/dt.
        output[3] = (1.0 - z_by)*lamba[BY] - mu*z_by; //z_by/dt.

        //dw/dt: Change to proportion of hosts immune to any strain sharing alleles with each strain.
        output[4] = (1.0 - w_ax)*(lamba[AX]+lamba[AY]+lamba[BX]) - mu*w_ax; //w_ax/dt.
        output[5] = (1.0 - w_ay)*(lamba[AY]+lamba[AX]+lamba[BY]) - mu*w_ay; //w_ay/dt.
        output[6] = (1.0 - w_bx)*(lamba[BX]+lamba[AX]+lamba[BY]) - mu*w_bx; //w_bx/dt.
        output[7] = (1.0 - w_by)*(lamba[BY]+lamba[AY]+lamba[BX]) - mu*w_by; //w_by/dt.

        //dy/dt: Change to proportions of hosts infectious with each strain.
        output[8] = ((1.0 - w_ax) + (1.0 - gamma)*(w_ax - z_ax))*lamba[AX] - sigma*y_ax; //z_ax/dt.
        output[9] = ((1.0 - w_ay) + (1.0 - gamma)*(w_ay - z_ay))*lamba[AY] - sigma*y_ay; //z_ay/dt.
        output[10] = ((1.0 - w_bx) + (1.0 - gamma)*(w_bx - z_bx))*lamba[BX] - sigma*y_bx; //z_bx/dt.
        output[11] = ((1.0 - w_by) + (1.0 - gamma)*(w_by - z_by))*lamba[BY] - sigma*y_by; //z_by/dt.
        */
    }

    //Returns the sum of force of infections for each strain overlapping strains[_strain] using the overlap matrix (including self).
    double phi(const unsigned int _strain, const std::vector<double>& _currentVals, const std::vector<double>& _betas) const
    {
        double output = 0.0;
        for (std::size_t i = 0; i<numStrains; i++)
        {
            if (overlapMatrix[_strain][i] >= 1) //If there is an overlap of at least 1, add the force of infection.
                output += _betas[i]*_currentVals[(numStrains*2)+i]; //FOI_i = beta*y_i.
        }

        return output;
    }

    //Calculate the overlap between each possible strain, storing the results in a matrix.
    void calculate_overlap_matrix()
    {
        //Generate list of all possible strains.
        std::vector<std::vector<unsigned int>> strains;
        recursive_permutations(numLoci, numAlleles, strains);

        //Calculate overlap matrix.
        for (std::size_t s=0; s<numStrains; s++) //For each strain s.
        {
            std::vector<unsigned int> row;
            for (std::size_t j=0; j<numStrains; j++) //Calculate number of overlaps for each other strain j.
            {
                //if (s == j) //Self-strain IS included so comment out.
                //{
                //    row.push_back(0)
                //    continue;
                //}

                //Calculate overlaps between strain s and strain j.
                unsigned int overlapCount = 0;
                for (std::size_t loci=0; loci<strains[s].size(); loci++) //Check for overlap at each loci.
                    if (strains[s][loci] == strains[j][loci])
                        overlapCount++;
                row.push_back(overlapCount);
            }
            overlapMatrix.push_back(row);
        }
    }

    void recursive_permutations(const unsigned int _n, const unsigned int _k, std::vector<std::vector<unsigned int>>& _strains, std::vector<unsigned int> _start = std::vector<unsigned int>(), const unsigned int _depth=0)
    {
        if (_depth == _n)
        {
            _strains.push_back(_start);
            return;
        }

        for (unsigned int i=0; i<_k; i++)
        {
            _start.push_back(i);
            recursive_permutations(_n, _k, _strains, _start, _depth+1);
            _start.pop_back();
        }
    }

    const unsigned int get_num_strains() const
    {
        return numStrains;
    }
};
