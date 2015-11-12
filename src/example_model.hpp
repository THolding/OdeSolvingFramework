#pragma once

class ExampleModel : public ModelDefinition
{
public:
    void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const
    {
        output[0] = -params[0]*currentValues[0]*currentValues[1] +params[1]*currentValues[1] + 1.0*params[2] - currentValues[0]*params[2];
        output[1] = +params[0]*currentValues[0]*currentValues[1] -params[1]*currentValues[1] - currentValues[1]*params[2];

        //std::cout << params[0] << "\t" << params[1] << "\t" << params[2] << "\t" << currentValues[0] << "\t" << currentValues[1] << "\n";
        //std::cout << "Total pop: " << currentValues[0]+currentValues[1] << "\n";
        //std::cout << "New infected: " << params[0]*currentValues[0]*currentValues[1] << "\n";
        //std::cout << "New recovered: " << +params[1]*currentValues[1] << "\n\n";
    }
};
