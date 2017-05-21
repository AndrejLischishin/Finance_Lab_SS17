//
//  integration_functions.cpp
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//
/** @file */
#include "integration_functions.hpp"
#include "random_functions.hpp"

double calculate_relative_error(double exact_value, double calculated_value)
{
	return fabs((exact_value-calculated_value)/exact_value);
}

/**
 * Simulates a wiener process.
 *
 * @param r Pointer to the gsl_rng object for generating standard normal distributed numbers
 * @param T Time period of simulated process
 * @param delta_t Step of discretisation
 *
 * @return Pointer to the vector of values at discretisation points
 */
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


void trap_rule(std::vector<double>* nodes, std::vector<double>* weights, int l)
{
	int Nl = pow(2,l)-1;
	double weight = (double) 3/(2*(Nl+1));
	double node;
    weights->push_back(weight);
   
	for(int i=1; i<=Nl; i++){
		node = (double) i/(Nl+1);
		weight = (double) 1/(Nl+1);
		nodes->push_back(node);
		weights->push_back(weight);
	}
	weight = (double) 3/(2*(Nl+1));
	weights->push_back(weight);
}

void clenshaw_curtis(std::vector<double>* nodes, std::vector<double>* weights, int l)
{
	int Nl = pow(2,l)-1;
	for(int i=1; i<Nl; i++)
	{
		nodes->push_back(.5*(1.-cos((double)(M_PI*i/(Nl+1.)))));
		double sum = 0.;

		for(int j=1; j<=(Nl+1)/2; j++)
		{
			sum = sum + (double) 1./(2.*j-1.)*sin((double)(2.*j-1.)*M_PI*i/(Nl+1.));
		}
		weights->push_back((double)2./(Nl+1.)*sin((double)M_PI*i/(Nl+1.))*sum);
	}
}

void gauss_legendre(std::vector<double>* nodes, std::vector<double>* weights, size_t l)
{	
	int Nl = pow(2,l)-1;

	double* wi = new double[1];
	double* xi = new double[1];
	
	double a=0.; 
	double b=1.;
	
	gsl_integration_glfixed_table* table;	
	table =	gsl_integration_glfixed_table_alloc(Nl);

	for(int j=0; j<Nl; j++)
	{
		gsl_integration_glfixed_point(a, b, j, xi, wi, table);
		nodes->push_back(xi[0]);
		weights->push_back(wi[0]);
	}

	gsl_integration_glfixed_table_free(table);
}

void monte_carlo(std::vector<double>* nodes, std::vector<double>* weights, int l, gsl_rng* r)
{		
	int Nl = pow(2,l)-1;
	for(int i=1; i<=Nl; i++)
	{
		nodes->push_back(random_number_01_GSL(r));
		weights->push_back(1./Nl);
	}
}

double call_option_integrand(double x, double s0, double mu, double sigma, double T, double K)
{
	double result = s0*exp((mu-0.5*sigma*sigma)*T+sigma*sqrt(T)*x)-K;
	if(result > 0.0)
		return result;
	else
		return 0.0;
}
