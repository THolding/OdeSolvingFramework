#pragma once

class GuptaMultistrainW : public ModelDefinition
{
public:
    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        const std::vector<double> beta(params.begin(), params.end()-3);
        const double gamma = params.at(params.size()-3);
        const double sigma = params.at(params.size()-2);
        const double mu = params.at(params.size()-1);

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
    }
};
