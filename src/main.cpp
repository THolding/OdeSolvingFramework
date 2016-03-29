#include <iostream>
#include <iomanip>
#include "model_driver.hpp"
#include "repertoire_multistrain_w.hpp"
#include "simulation_sets.hpp"

int main()
{
    /*const unsigned int numLoci = 2;
    const unsigned int numAlleles = 4;
    const double beta = 0.6;
    const double gamma = 0.70;
    const double sigma = 0.3;
    const double mu = 0.001;
    const double initInfectiousTotal = 0.001;

    run_diversity_sweep(numLoci, {numAlleles}, beta, gamma, sigma, mu, initInfectiousTotal);*/

    ///Num alleles sweep (diversity_sweep).
    /*const unsigned int numLoci = 2;
    const std::vector<unsigned int> numAlleles = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    const double beta = 0.6;
    const double gamma = 0.70;
    const double sigma = 0.3;
    const double mu = 0.001;
    const double initInfectiousTotal = 0.001;
    run_diversity_sweep(numLoci, numAlleles, beta, gamma, sigma, mu, initInfectiousTotal, "diversity_sweep_test");*/

    ///Randomly selected strains, num strains sweep.
    const unsigned int numLoci = 2;
    const unsigned int numAlleles = 11; //105 possible strains.
    const std::vector<unsigned int> numStrainsList = { 2, 4, 8, 12, 16, 20, 25, 30, 40, 55};
    const unsigned int numReps = 30;
    const double beta = 0.6;
    const double gamma = 0.70;
    const double sigma = 0.3;
    const double mu = 0.001;
    const double initInfectiousTotal = 0.001;
    run_random_initial_strain_sweep(numLoci, numAlleles, numStrainsList, numReps, beta, gamma, sigma, mu, initInfectiousTotal, "random_num_strains_sweep");


    /*RepertoireMultistrainW modelDef(2,4);

    const double beta = 0.9;

    std::vector<double> params;
    params.push_back(0.6); //beta_1
    params.push_back(0.6); //beta_2
    params.push_back(0.6); //beta_3
    params.push_back(0.6); //beta_4
    params.push_back(0.6); //beta_5
    params.push_back(0.6); //beta_6
    params.push_back(0.8); //gamma (cross-reactivity)
    params.push_back(0.3); //sigma
    params.push_back(0.001); //mu
    */

    /*for (std::size_t i=0; i<modelDef.get_num_strains(); i++)
        params.push_back(0.9); //beta_i.
    params.push_back(0.6); //gamma (cross-immunity).
    params.push_back(0.3); //sigma (recovery).
    params.push_back(0.0001); //mu (births and deaths).*/

    /*for (std::size_t i=0; i<modelDef.get_num_strains(); i++)
        params.push_back(0.2); //beta_i.
    params.push_back(0.6); //gamma (cross-immunity).
    params.push_back(0.1); //sigma (recovery).
    params.push_back(0.001); //mu (births and deaths).*/

    /*std::vector<double> ys = { 0.01, 0.01, 0.01, 0.00999, 0.01, 0.01001 };
    std::vector<double> init;
    modelDef.calculate_initial_values(ys, params, init);

    //for (auto i:init)
    //    std::cout << i << "\n";

    ModelDriver model(&modelDef, params, init);
    model.set_dt(0.05);
    model.set_max_time(50000);
    model.set_output_frequency(20);
    //model.set_stop_threshold(0);
    model.run("repertoire_multistrain_w");
    modelDef.export_num_strains("repertoire_multistrain_w");

    model.export_output();*/

    return 0;
}
