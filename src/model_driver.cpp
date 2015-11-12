#include "model_driver.hpp"
#include <cmath>
#include "runge_kutta_ck_solver.hpp"
#include "utilities.hpp"

//Solves the system of ODEs.
void ModelDriver::run(std::string _runName)
{
    runName = _runName;
    stopCondition = StopCondition(stopThreshold, maxTime);
    RungeKuttaCKSolver solver(this);
    solver.solve();
}

//Updates the model, i.e. commits proposed change.
void ModelDriver::update(const double _t, const std::vector<double>& _currentValues)
{
    curVals = _currentValues; //Update current value.

    //Log to output file.
    std::vector<double> newRow;
    newRow.push_back(_t);
    newRow.insert(newRow.end(), _currentValues.begin(), _currentValues.end());
    output.push_back(newRow);
}

//Returns a copy of the output.
std::vector<std::vector<double>> ModelDriver::get_output() const
{
    return output;
}

//Export output to file.
void ModelDriver::export_output() const
{
    matrixToFile(output, runName+".csv", ", ");
}


//Internal logic of stop condition check. Returns true when stop condition is met.
bool StopCondition::check(const double _time, const std::vector<double>& _state)
{
    if (_time >= maxTime)
        return true;

    if (!first) //If not the first check.
    {
        const double dt = _time - lastTime;

        //If any one compartment changes more than threshold then it is not time to stop.
        for (unsigned int i=0; i<_state.size(); i++)
        {
            double diff = std::abs(lastState[i] - _state[i]);
            if (diff > threshold*dt)
            {
                lastTime = _time; //Update last time and state before returning.
                lastState = _state;
                return false;
            }
        }

        //If program gets to here then no compartments have changed enough and stop condition is true.
        return true;
    }
    else //First check so no previous values.
    {
        lastTime = _time;
        lastState = _state;
        first = false;
        return false;
    }
}
