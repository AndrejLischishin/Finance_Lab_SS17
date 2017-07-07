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

#endif /* exotic_options.hpp */
