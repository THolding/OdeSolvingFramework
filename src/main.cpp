#include <iostream>
#include "model_driver.hpp"
#include "gupta_multistrain_w.hpp"

int main()
{
    GuptaMultistrainW modelDef(4,3);

    std::vector<double> params;
    for (std::size_t i=0; i<modelDef.get_num_strains(); i++)
        params.push_back(0.5); //beta_i.
    params.push_back(0.55); //gamma (cross-immunity).
    params.push_back(0.2); //sigma (recovery).
    params.push_back(0.0008); //mu (births and deaths).

    std::vector<double> init;
    for (std::size_t i=0; i<modelDef.get_num_strains(); i++)
        init.push_back(0.0); //z_i.

    for (std::size_t i=0; i<modelDef.get_num_strains(); i++)
        init.push_back(0.0); //w_i.

    for (std::size_t i=0; i<modelDef.get_num_strains()-1; i++)
        init.push_back(0.001); //y_i.
    init.push_back(0.001001); //y_(n-1).

    /*ModelDriver model(&modelDef, params, init);
    model.set_dt(0.05);
    model.set_max_time(10000);
    model.set_output_frequency(10);
    model.run("gupta_multistrain_w");

    model.export_output();*/

    return 0;
}
