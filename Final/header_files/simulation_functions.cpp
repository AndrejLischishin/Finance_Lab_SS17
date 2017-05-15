//
//  simulation_functions.cpp
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//
/** @file */
#include "simulation_functions.hpp"

/**
 * Simulates a wiener process.
 *
 * @param r Pointer to the gsl_rng object for generating standard normal distributed numbers
 * @param T Time period of simulated process
 * @param delta_t Step of discretisation
 *
 * @return Pointer to the vector of values at discretisation points
 */
std::vector<double>* wiener_process(gsl_rng* r, double T, double delta_t)
{
    int M = (int)(T/delta_t);
    std::vector<double> *w = new std::vector<double>(M+1);
    (*w)[0] = 0.0;

    for(int i=0; i<M; i++)
    {
        (*w)[i+1] = (*w)[i]+sqrt(delta_t)*gsl_ran_ugaussian(r);
    }
    
    return w;
}

/**
 * Simulates brownian motion path for the given values of wiener process.
 *
 * @param r Pointer to the gsl_rng object for generating standard normal distributed numbers
 * @param T Time period of simulated process
 * @param delta_t Step of discretisation
 * @param w Pointer to the vector with values of wiener process at discretisation points
 * @param s0 Value of brownian_motion at time = 0
 * @param mu Drift
 * @param sigma Volatility
 *
 * @return Pointer to the vector of values at discretisation points
 */
std::vector<double>* brownian_motion(gsl_rng* r, double T, double delta_t, std::vector<double>* w, double s0, double mu, double sigma)
{
    int M = (int)(T/delta_t);
    std::vector<double> *s = new std::vector<double>(M+1);
    (*s)[0] = s0;
    for(int i=0; i<M; i++)
    {
        (*s)[i+1] = s0*exp((mu-0.5*pow(sigma, 2.0))*(i+1)*delta_t+sigma*(*w)[i+1]);
    }
    
    return s;
}
