#include <iostream>
#include <iomanip>
#include "model_driver.hpp"
#include "example_model.hpp"

int main()
{
    ExampleModel model;

    std::vector<double> params =
    {
        1000, //H = 1000, //host population size
        10000, //M = 10000, //mosquito population size
        0.33, //a = 0.33, //Mosquito bite rate to humans
        0.5, //b = 0.5, //probability infectious bite leads to infection in humans
        0.25, //c = 0.25, //probability that a mosquito biting infectious human becomes infectious
        0.15, //mu_m = 0.15, //Mosquito death rate (1-p)
        0.0001, //mu_h = 0.0001, //host death rate
        0.2, //tau = 0.2, //rate infected mosquitoes become infectious
        0.02 //sigma = 0.02 //recovery rate of hosts.
    };

    std::vector<double> init =
    {
        0.01, //Z
        0.00, //X
        0.01 //Y
    };

    ModelDriver driver(&model, params, init);
    driver.set_output_frequency(100);
    driver.run("test");
    driver.export_output();

    return 0;
}
