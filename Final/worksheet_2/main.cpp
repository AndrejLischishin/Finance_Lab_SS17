//
//  main.cpp
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//
/** @file */
#include "../header_files/random_functions.hpp"
#include "../header_files/simulation_functions.hpp"

//using namespaces in order to be able to use same names more than once
namespace Task_1 {

  double s0;
  double mu;
  std::vector<double> sigma(5);
  std::vector<double> V_mean(5);
  std::vector<double>* w;
  std::vector<double>* s;
  double T;
  double delta_t;
  double K;


}

namespace Task_2 {

  double s0;
  double mu;
  double sigma;
  std::vector<double> delta_t(4);
  std::vector<double> V_mean(4);
  std::vector<double> V_variance(4);
  std::vector<double>* w;
  std::vector<double>* s;
  double T;
  double K;

}

/**
 * Main function to run all exercises of worksheet 2.
 *
 * @param argc Integer argument for main
 * @param argv Char array argument for main
 *
 * @return Returns 0 if everything worked fine
 *
 */
int main(int argc, char* argv[]){

	//std::cout << "Prepared everything for worksheet 2." << std::endl;

  gsl_rng* r;

  //seeding
  unsigned long seed = time(NULL);
  //memory allocation
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, seed);


    //////////////////////////////////////////////////////////
    //////////////////////////Task_1//////////////////////////
    //////////////////////////////////////////////////////////

    Task_1::s0 = 10;
    Task_1::mu = 0.1;
    Task_1::sigma = {0.0, 0.2, 0.4, 0.6, 0.8};  // different sigmas
    Task_1::T = 2;                                // time intervall T
    Task_1::delta_t = 0.2;
    Task_1::K = 10;                               // strike price of the option

    double N = 1000;                             // number of simulations


    for (int i = 0; i < Task_1::sigma.size(); i++) {
      for (int j = 0; j < N; j++) {
        Task_1::w = wiener_process(r,Task_1::T, Task_1::delta_t);
        Task_1::s = brownian_motion(r,Task_1::T, Task_1::delta_t,
           Task_1::w, Task_1::s0, Task_1::mu, Task_1::sigma[i]);
        Task_1::V_mean[i] += std::max((*Task_1::s).back()-Task_1::K,0.0);
      }
      Task_1::V_mean[i] = Task_1::V_mean[i]/1000.0;
    }

//Plot V aginst sigma
//Want to use Python for this, but there are still some linker problems
//Working no it

    //////////////////////////////////////////////////////////
    //////////////////////////Task_2//////////////////////////
    //////////////////////////////////////////////////////////

    Task_2::s0 = 10;
    Task_2::mu = 0.1;
    Task_2::T = 2;                                // time intervall T
    Task_2::delta_t = {0.2,0.8,1.,2.};
    Task_2::K = 10;                               // strike price of the option

    std::vector<double>* helper_1 = new std::vector<double>;
    double helper_2;

    for (int i = 0; i < Task_2::delta_t.size(); i++) {
      for (int j = 0; j < N; j++) {
        Task_2::w = wiener_process(r,Task_2::T, Task_2::delta_t[i]);
        Task_2::s = brownian_motion(r,Task_2::T, Task_2::delta_t[i],
           Task_2::w, Task_2::s0, Task_2::mu, Task_2::sigma);
        helper_2 = std::max((*Task_2::s).back()-Task_2::K,0.0);
        Task_2::V_mean[i] += helper_2;
        helper_1->push_back(helper_2);
      }
      Task_2::V_mean[i] = Task_2::V_mean[i]/1000.0;
      Task_2::V_variance[i] = sigma_algorithm(helper_1,N);
    }

//Plot V_variance against delta_t
//Want to use Python for this, but there are still some linker problems
//Working on it

    gsl_rng_free(r);
    return 0;
}
