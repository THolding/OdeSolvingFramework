#pragma once
#include <string>
#include <vector>
#include "model_definition.hpp"

class StopCondition
{
private:
    double threshold;
    double maxTime;
    std::vector<double> lastState;
    double lastTime;
    bool first=true;
public:
    StopCondition(const double _threshold = 1.0, const double _maxTime = 1.0) : threshold(_threshold), maxTime(_maxTime) {  }
    bool check(const double _time, const std::vector<double>& _state);
};

class ModelDriver
{
private:
    ModelDefinition* modelDefinition;
    StopCondition stopCondition;
    std::vector<double> params;
    std::vector<double> curVals;

    std::string runName;
    double startTime;
    double maxTime;
    double dt;
    double stopThreshold;

    std::vector<std::vector<double>> output;
public:
    ModelDriver(ModelDefinition* _modelDefinition, std::vector<double> _params, std::vector<double> _initVals)
        : modelDefinition(_modelDefinition), params(_params), curVals(_initVals),
          runName("default_run"), startTime(0.0), maxTime(100.0), dt(0.0005), stopThreshold(0.000001)
    {  }
    void run(std::string _runName = "default_runName"); //Solves the system of ODEs.
    void update(const double _t, const std::vector<double>& _currentValues); //Updates the model, i.e. commits proposed change.
    bool stop_condition_met(const double _time, const std::vector<double>& _state) { return stopCondition.check(_time, _state); }; //Returns true if stop condition has been met.
    std::vector<std::vector<double>> get_output() const; //Returns a copy of the output.
    void export_output() const; //Exports to csv file.

    //Model settings / getter-setters.
    void set_start_time(const double _startTime) { startTime = _startTime; }
    void set_max_time(const double _maxTime) { maxTime = _maxTime; }
    void set_dt(const double _dt) { dt = _dt; }
    void set_stop_threshold(const double _stopThreshold) { stopThreshold = _stopThreshold; }
    double get_start_time() const { return startTime; }
    double get_max_time() const { return maxTime; }
    double get_dt() const { return dt; }
    double get_num_compartments() const { return curVals.size(); }
    const ModelDefinition* get_model_definition() const { return modelDefinition; }
    std::vector<double> get_params() const { return params; }
    std::vector<double> get_current_values() const { return curVals; }
    //
};
