//
//  main.cpp
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//
/** @file */
#include "../header_files/random_functions.hpp"
#include "../header_files/simulation_functions.hpp"
#include "../header_files/integration_functions.hpp"

template<typename... Args>
double integrate_by_point_evaluation(double (*function_to_integrate)(double x, Args... rest), std::vector<double>* nodes, std::vector<double>* weights, int n, Args... rest)
{
	double result = 0.0;
	for(int i=0; i<n; i++)
	{
		result += (*weights)[i]*function_to_integrate((*nodes)[i], rest...);
	}

	return result;
}

double function_to_integrate(double x)
{
	return x;
}

double f_gamma(double x, double gamma)
{
	return 1.0+gamma*exp(0.5*x);
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
    
	std::cout << "Prepared everything for worksheet 2." << std::endl;

	int l = 10;
	
	gsl_rng* r;
    
    //seeding
    unsigned long seed = time(NULL);
    //memory allocation
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);

	std::vector<double>* nodes;
    nodes = new std::vector<double>;
	std::vector<double>* weights;
    weights = new std::vector<double>;
	trap_rule(nodes, weights, l);
	std::cout << "Trapezoidal rule:" << std::endl;
	std::cout << integrate_by_point_evaluation(function_to_integrate, nodes, weights, pow(2,l)-1) << std::endl;
	std::cout << integrate_by_point_evaluation(f_gamma, nodes, weights, pow(2,l)-1, 1.) << std::endl;
	std::cout << integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, 10., 0.1, 0.2, 1., 10.) << std::endl;
	nodes->clear();
	weights->clear();
	clenshaw_curtis(nodes, weights, l);
	std::cout << "Clenshaw curtis:" << std::endl;
	std::cout << integrate_by_point_evaluation(function_to_integrate, nodes, weights, pow(2,l)-1) << std::endl;
	std::cout << integrate_by_point_evaluation(f_gamma, nodes, weights, pow(2,l)-1, 1.) << std::endl;
	std::cout << integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, 10., 0.1, 0.2, 1., 10.) << std::endl;
	nodes->clear();
	weights->clear();
	monte_carlo(nodes, weights, pow(2,l)-1, r);
	std::cout << "Monte Carlo:" << std::endl;
	std::cout << integrate_by_point_evaluation(function_to_integrate, nodes, weights, pow(2,l)-1) << std::endl;
	std::cout << integrate_by_point_evaluation(f_gamma, nodes, weights, pow(2,l)-1, 1.) << std::endl;
	std::cout << integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, 10., 0.1, 0.2, 1., 10.) << std::endl;
	nodes->clear();
	weights->clear();
	gauss_legendre(nodes, weights, pow(2,l)-1);
	std::cout << "Gauss Legendre:" << std::endl;
	std::cout << integrate_by_point_evaluation(function_to_integrate, nodes, weights, pow(2,l)-1) << std::endl;
	std::cout << integrate_by_point_evaluation(f_gamma, nodes, weights, pow(2,l)-1, 1.) << std::endl;
	std::cout << integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, 10., 0.1, 0.2, 1., 10.) << std::endl;
	
    return 0;
}
