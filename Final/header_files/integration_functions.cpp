//
//  integration_functions.cpp
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//
/** @file */
#include "integration_functions.hpp"
#include "random_functions.hpp"

/**
 * Returns the relative error between the given exact value and a calculated value.
 *
 * @param exact_value The exact value of a calculation
 * @param calculated_value The calculated value of an algorithm
 *
 * @return The relative error of the algorithm
 */
double calculate_relative_error(double exact_value, double calculated_value)
{
	return fabs((exact_value-calculated_value)/exact_value);
}


/**
 * Evaluates a function at given points and multiplies every value with a given weight and sum all the results up.
 *
 * @param (*function_to_integrate)(double x, Args... rest) The function with its additional parameters the will be integrated
 * @param nodes The points on which the function will be evaluated
 * @param weights The weights which will be multiplied with the values of the function at the different points
 * @param n The number of evaluation-points
 * @param Args... The additional arguments for the function
 *
 * @return The value of the sum (which approximates the integral of the function on a certain interval)
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


/**
 * Calculates nodes and weights by the trapezoidal rule on the interval \f$[0,1]\f$.
 *
 * @param nodes The vector that will include the points
 * @param weights The vector that will include all the weights at the different points
 * @param l The number of points will be calculated by \f$2^l-1\f$
 */
void trap_rule(std::vector<double>* nodes, std::vector<double>* weights, int l)
{
	unsigned int Nl = pow(2,l)-1;
	double weight = 3./(double)(2*(Nl+1));
	double node;
    weights->push_back(weight);

	for(int i=1; i<=Nl; i++){
		node = (double)i/(double)(Nl+1);
		nodes->push_back(node);
		if (i>=2&&i<Nl) {
			weight = 1./(double)(Nl+1);
			weights->push_back(weight);
		}
	}
	weight = 3./(double)(2*(Nl+1));
	weights->push_back(weight);
}


/**
 * Calculates nodes and weights by the Clenshaw-Curtis-Rule on the interval \f$[0,1]\f$.
 *
 * @param nodes The vector that will include the points
 * @param weights The vector that will include all the weights at the different points
 * @param l The number of points will be calculated by \f$2^l-1\f$
 */
void clenshaw_curtis(std::vector<double>* nodes, std::vector<double>* weights, int l)
{
	unsigned int Nl = pow(2,l)-1;

	for(unsigned int i=1; i<=Nl; i++)
	{
		nodes->push_back(.5*(1.-cos((double)(M_PI*i/(Nl+1.)))));
		double sum = 0.;

		for(unsigned int j=1; j<=(Nl+1)/2; j++)
		{
			sum = sum + (double) (1./(2.*j-1.)) * sin( (double)(2.*j-1.) * M_PI * (i/(double)(Nl+1.)) );
		}
		weights->push_back( (double)(2./(Nl+1.)) * sin( M_PI * ((double)i/(double)(Nl+1.)) ) * sum );
	}


}


/**
 * Calculates nodes and weights by the Gauss-Legendre-Rule on the interval \f$[0,1]\f$.
 *
 * @param nodes The vector that will include the points
 * @param weights The vector that will include all the weights at the different points
 * @param l The number of points will be calculated by \f$2^l-1\f$
 */
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


/**
 * Calculates nodes and weights by the Monte-Carlo-Approach on the interval \f$[0,1]\f$.
 *
 * @param nodes The vector that will include the points
 * @param weights The vector that will include all the weights at the different points
 * @param l The number of points will be calculated by \f$2^l-1\f$
 */
void monte_carlo(std::vector<double>* nodes, std::vector<double>* weights, int l, gsl_rng* r)
{
	int Nl = pow(2,l)-1;
	for(int i=1; i<=Nl; i++)
	{
		nodes->push_back(random_number_01_GSL(r));
		weights->push_back(1./Nl);
	}
}


/**
 * Calculates the value of the integrand for european call options.
 *
 * @param x Point to evaluate the integrand at
 * @param s0 Value of \f$S\f$ at point \f$0\f$
 * @param mu Value of \f$\mu\f$
 * @param sigma Value of \f$\sigma\f$
 * @param T Right bound of the time interval
 * @param K Strike price
 *
 * @return The calculated value
 */
double call_option_integrand(double x, double s0, double mu, double sigma, double T, double K)
{
	x = normal_inverse_cdf(x);
	double result = (s0*exp(((mu-0.5*sigma*sigma)*T)+(sigma*sqrt(T)*x)))-K;
	if(result > 0.0)
		return result;
	else
		return 0.0;
}
