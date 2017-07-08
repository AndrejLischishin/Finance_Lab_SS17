//
//  exotic_options.hpp
//

#ifndef exotic_options_hpp
#define exotic_options_hpp

#include <vector>
#include <math.h>
#include "random_functions.hpp"

double payoff_discrete_down_out_call(std::vector<double>* x, double s0, double r, double T, int M, double K, double sigma, double B);
double payoff_discrete_lookback(std::vector<double>* x, double s0, double r, double T, int M, double K, double sigma);
double black_scholes_down_out_call(double s0, double K, double T, double sigma, double r, double B);
double d_S_K(double S, double K, double r, double sigma, double T);
double V_bs(double S, double K, double r, double sigma, double T);

#endif /* exotic_options.hpp */
