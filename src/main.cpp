#include <iostream>
#include "model_driver.hpp"
#include "repertoire_multistrain_w.hpp"

int main()
{
    RepertoireMultistrainW modelDef(2,4);

    std::vector<double> params;
    for (std::size_t i=0; i<modelDef.get_num_strains(); i++)
        params.push_back(0.9); //beta_i.
    params.push_back(0.6); //gamma (cross-immunity).
    params.push_back(0.3); //sigma (recovery).
    params.push_back(0.0001); //mu (births and deaths).
    /*for (std::size_t i=0; i<modelDef.get_num_strains(); i++)
        params.push_back(0.2); //beta_i.
    params.push_back(0.6); //gamma (cross-immunity).
    params.push_back(0.1); //sigma (recovery).
    params.push_back(0.001); //mu (births and deaths).*/

    std::vector<double> init;
    for (std::size_t i=0; i<modelDef.get_num_strains()-1; i++)
        init.push_back(0.001); //z_i.
    init.push_back(0.001001); //z_(n-1).

    const double w = 0.005;
    for (std::size_t i=0; i<modelDef.get_num_strains(); i++)
        init.push_back(w); //w_i.

    for (std::size_t i=0; i<modelDef.get_num_strains()-1; i++)
        init.push_back(0.001); //y_i.
    init.push_back(0.001001); //y_(n-1).

    std::vector<double> temp = init;
    modelDef.calc_derivatives(init, temp, params);
    for (auto i:temp)
        std::cout << i << ", ";
    std::cout << "\n\n";



    /*ModelDriver model(&modelDef, params, init);
    model.set_dt(0.05);
    model.set_max_time(1);
    model.set_output_frequency(1);
    model.set_stop_threshold(0);
    model.run("repertoire_multistrain_w");
    modelDef.export_num_strains("repertoire_multistrain_w");

    model.export_output();*/

    return 0;
}
