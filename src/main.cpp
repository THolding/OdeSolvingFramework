#include <iostream>
#include "model_driver.hpp"
#include "gupta_multistrain_w.hpp"

int main()
{
    std::vector<double> params;
    params.push_back(0.5); //beta_ax.
    params.push_back(0.5); //beta_ax_ay.
    params.push_back(0.5); //beta_ax_bx.
    params.push_back(0.5); //beta_ax_by.
    params.push_back(0.55); //gamma (cross-immunity).
    params.push_back(0.2); //sigma (recovery).
    params.push_back(0.0008); //mu (births and deaths).

    std::vector<double> init;
    init.push_back(0.0); //z_ax.
    init.push_back(0.0); //z_ay.
    init.push_back(0.0); //z_bx.
    init.push_back(0.0); //z_by.
    init.push_back(0.0); //w_ax.
    init.push_back(0.0); //w_ay.
    init.push_back(0.0); //w_bx.
    init.push_back(0.0); //w_by.
    init.push_back(0.001); //y_ax.
    init.push_back(0.001); //y_ay.
    init.push_back(0.001); //y_bx.
    init.push_back(0.001001); //y_by.

    GuptaMultistrainW modelDef;

    ModelDriver model(&modelDef, params, init);
    model.set_dt(0.05);
    model.set_max_time(10000);
    model.set_output_frequency(10);
    model.run("gupta_multistrain_w");

    model.export_output();

    return 0;
}
