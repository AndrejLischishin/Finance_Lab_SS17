//
//  integration_functions.hpp
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//

#ifndef integration_functions_hpp
#define integration_functions_hpp

#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>

double calculate_relative_error(double exact_value, double calculated_value);
template<typename... Args>
double integrate_by_point_evaluation(double (*function_to_integrate)(double x, Args... rest), std::vector<double>* points, std::vector<double>* weights, int n, Args... rest);
void trap_rule(std::vector<double>* nodes, std::vector<double>* weights, int l);
void clenshaw_curtis(std::vector<double>* nodes, std::vector<double>* weights, int l);
void gauss_legendre(std::vector<double>* nodes, std::vector<double>* weights, size_t l);
void monte_carlo(std::vector<double>* nodes, std::vector<double>* weights, int l, gsl_rng* r);
double call_option_integrand(double x, double s0, double mu, double sigma, double T, double K);

#endif /* integration_functions_hpp */
