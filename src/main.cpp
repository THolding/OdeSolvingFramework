#include <iostream>
#include "model_driver.hpp"
#include "example_model.hpp"

int main()
{
    ExampleModel modelDef;

    ModelDriver model(&modelDef, {2.0, 1.0, 0.0}, {0.99, 0.01});
    model.set_dt(0.001);
    model.set_max_time(100);
    model.run();

    model.export_output();

    return 0;
}
