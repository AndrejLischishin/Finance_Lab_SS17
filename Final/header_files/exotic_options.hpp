//
//  exotic_options.hpp
//

#ifndef exotic_options_hpp
#define exotic_options_hpp

#include <vector>
#include <math.h>
#include "random_functions.hpp"

double payoff_discrete_down_out_call(std::vector<double> x, double s0, double r, double T, int M, double K, double sigma, double B);
double payoff_discrete_lookback(std::vector<double> x, double s0, double r, double T, int M, double K, double sigma);
double asian_option_call_integrand_arithmetic(std::vector<double>* x,double S0,double K, double sigma, double mu,int M, double T, bool use_bb);
double asian_option_call_integrand_arithmetic_control_variates(std::vector<double>* x,double S0,double K, double sigma, double mu,int M, double T, bool use_bb);
double lookback_call_integrand_fixed(std::vector<double>* x,double S0,double K, double sigma, double mu,int M, double T, bool use_bb);
double barrier_integrand(std::vector<double>* x,double S0,double K, double sigma, double mu,int M, double T, double B, bool use_bb);
#endif /* exotic_options.hpp */
