#pragma once
#include <vector>

//Interface allowing polymorphism with different model definitions.
class ModelDefinition
{
public:
    //Performs one preliminary step (called multiple times before 'committing' a step).
    virtual void calc_derivatives(std::vector<double>& currentValues, std::vector<double>& output, const std::vector<double>& params) const=0;
};
