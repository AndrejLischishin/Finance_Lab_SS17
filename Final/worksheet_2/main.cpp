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
//plot with Python in C++ programs
//#include "../header_files/matplotlibcpp.hpp"



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

double call_option_exact_expected_value(double s0, double mu, double T, double sigma, double K)
{
	double chi = (log(K/s0)-((mu-sigma*sigma*0.5)*T))/(sigma*sqrt(T));
	return s0*exp(mu*T)*normal_cdf(sigma*sqrt(T)-chi)-K*normal_cdf(-chi);
}

//using namespaces in order to be able to use same names more than once
namespace Task_1 
{
  double s0;
  double mu;
  std::vector<double> sigma(5);
  std::vector<double> V_mean(5);
  std::vector<double>* w;
  std::vector<double>* s;
  double T;
  double delta_t;
  double K;
  double N;
}

namespace Task_2 
{
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
  double N;
}

namespace Task_10
{
  double s0;
  double mu;
  double sigma;
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
 */

//Plot, but because with Python and system depending, comented
//namespace plt = matplotlibcpp;
int main(int argc, char* argv[])
{
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
  	Task_1::T = 2;                              // time intervall T
  	Task_1::delta_t = 0.2;
  	Task_1::K = 10;                             // strike price of the option

  	Task_1::N = 1000.; 										 			// number of simulations

	//allocation of memmory
	Task_1::w = new std::vector<double>;
	Task_1::s = new std::vector<double>;


  	for(unsigned int i = 0; i < Task_1::sigma.size(); i++) {
    	for(int j = 0; j < Task_1::N; j++) {

			//simulating wiener_process
      		Task_1::w = wiener_process(r,Task_1::T, Task_1::delta_t);
			//simulating brownian_motion
			Task_1::s = brownian_motion(r,Task_1::T, Task_1::delta_t,
    		Task_1::w, Task_1::s0, Task_1::mu, Task_1::sigma[i]);
			//calculating of the Payoff on the maturity day
			Task_1::V_mean[i] += std::max((*Task_1::s).back()-Task_1::K,0.0);
    	}

		//calculating mean for different sigma
    	Task_1::V_mean[i] = Task_1::V_mean[i]/Task_1::N;

		//clearing vectors for next iteration
		Task_1::w->clear();
		Task_1::s->clear();
	}

//Plot V aginst sigma, but because with Python and system depending, comented
//plt::plot( Task_1::sigma, Task_1::V_mean, "r-");
//plt::save("./Task_1.png");
//plt::show();


	//////////////////////////////////////////////////////////
    //////////////////////////Task_2//////////////////////////
    //////////////////////////////////////////////////////////

    Task_2::s0 = 10;
    Task_2::mu = 0.1;
	Task_2::sigma = 0.2;
    Task_2::T = 2;                      // time intervall T
    Task_2::delta_t = {0.2,0.4,1.,2.};
    Task_2::K = 10.;											// strike price of the option

    Task_2::N = 1000.;									// number od simulations


	// helpvector allocation of the memmory
    std::vector<double>* helper_1 = new std::vector<double>;
	// helpobject for some technicall further calculations
	double helper_2 = 0;
	//allocation of memmory
	Task_2::w = new std::vector<double>;
	Task_2::s = new std::vector<double>;


	for(unsigned int i = 0; i < Task_2::delta_t.size(); i++) {
	    for(int j = 0; j < Task_2::N; j++) {

			// simulating wiener_process
			Task_2::w = wiener_process(r,Task_2::T, Task_2::delta_t[i]);
			//simulating brownian_motion
    	    Task_2::s = brownian_motion(r,Task_2::T, Task_2::delta_t[i], Task_2::w, Task_2::s0, Task_2::mu, Task_2::sigma);
			//calculating of the Payoff on the maturity day
    	    helper_2 = std::max((*Task_2::s).back()-Task_2::K,0.0);
			//helper for the calculating of the mean
    	    Task_2::V_mean[i] += helper_2;
			//vector v of Payoffs for fixed delta_t, simulated 1000times=>v.length=1000
    	    helper_1->push_back(helper_2);
		}

		//calculating mean for different delta_t
		Task_2::V_mean[i] = Task_2::V_mean[i]/Task_2::N;
		// calculating sigma(variance) for different delta_t
		Task_2::V_variance[i] = sigma_algorithm(helper_1,Task_2::N);

		// clearing vectors for next iteration
		Task_2::w->clear();
		Task_2::s->clear();
		helper_1->clear();
	}


//Plot V_variance against delta_t, same reason as above
//plt::plot(Task_2::delta_t, Task_2::V_variance, "r-");
//plt::save("./Task_2.png");
//plt::show();




	int l = 10;
	int l_max = 10;	

	std::vector<double>* nodes;
    nodes = new std::vector<double>;
	std::vector<double>* weights;
    weights = new std::vector<double>;

/*
	trap_rule(nodes, weights, l);
	std::cout << "Trapezoidal rule:" << std::endl;
	std::cout << integrate_by_point_evaluation(function_to_integrate, nodes, weights, l) << std::endl;
	std::cout << integrate_by_point_evaluation(f_gamma, nodes, weights, l, 1.) << std::endl;
	std::cout << integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, 10., 0.1, 0.2, 1., 10.) << std::endl;
	nodes->clear();
	weights->clear();
	clenshaw_curtis(nodes, weights, l);
	std::cout << "Clenshaw curtis:" << std::endl;
	std::cout << integrate_by_point_evaluation(function_to_integrate, nodes, weights, l) << std::endl;
	std::cout << integrate_by_point_evaluation(f_gamma, nodes, weights, l, 1.) << std::endl;
	std::cout << integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, 10., 0.1, 0.2, 1., 10.) << std::endl;
	nodes->clear();
	weights->clear();
	monte_carlo(nodes, weights, l, r);
	std::cout << "Monte Carlo:" << std::endl;
	std::cout << integrate_by_point_evaluation(function_to_integrate, nodes, weights, l) << std::endl;
	std::cout << integrate_by_point_evaluation(f_gamma, nodes, weights, l, 1.) << std::endl;
	std::cout << integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, 10., 0.1, 0.2, 1., 10.) << std::endl;
	nodes->clear();
	weights->clear();
	gauss_legendre(nodes, weights, l);
	std::cout << "Gauss Legendre:" << std::endl;
	std::cout << integrate_by_point_evaluation(function_to_integrate, nodes, weights, l) << std::endl;
	std::cout << integrate_by_point_evaluation(f_gamma, nodes, weights, l, 1.) << std::endl;
	std::cout << integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, 10., 0.1, 0.2, 1., 10.) << std::endl;
*/

	
	/* Integration of f_gamma by different formulas */
	std::ofstream myfile;
	myfile.open("output/relative_errors_f_gamma.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	double exact_result = 2.*exp(0.5)-1.;
	double calculated_result;
	for(l=1; l<=l_max; l++)
	{
		myfile<<pow(2,l)-1<<"	";

		nodes->clear();
		weights->clear();
		monte_carlo(nodes, weights, l, r);
		calculated_result = integrate_by_point_evaluation(f_gamma, nodes, weights, pow(2,l)-1, 1.);
		myfile<<calculate_relative_error(exact_result, calculated_result)<<"	";

		nodes->clear();
		weights->clear();
		trap_rule(nodes, weights, l);
		calculated_result = integrate_by_point_evaluation(f_gamma, nodes, weights, pow(2,l)-1, 1.);
		myfile<<calculate_relative_error(exact_result, calculated_result)<<"	";

		nodes->clear();
		weights->clear();
		clenshaw_curtis(nodes, weights, l);
		calculated_result = integrate_by_point_evaluation(f_gamma, nodes, weights, pow(2,l)-1, 1.);
		myfile<<calculate_relative_error(exact_result, calculated_result)<<"	";

		nodes->clear();
		weights->clear();
		gauss_legendre(nodes, weights, l);
		calculated_result = integrate_by_point_evaluation(f_gamma, nodes, weights, pow(2,l)-1, 1.);
		myfile<<calculate_relative_error(exact_result, calculated_result)<<std::endl;
	}
    myfile.close();


	/* Integration of the call option integrand by different formulas */
	myfile.open("output/relative_errors_K_10.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	Task_10::s0 = 10.;
	Task_10::mu = 0.1;
	Task_10::T = 1.;
	Task_10::sigma = 0.2;
	Task_10::K = 10.;
	exact_result = call_option_exact_expected_value(Task_10::s0, Task_10::mu, Task_10::T, Task_10::sigma, Task_10::K);
	std::cout << "Exact result: " << calculated_result << std::endl;
	for(l=1; l<=l_max; l++)
	{
		myfile<<pow(2,l)-1<<"	";

		/* Monte Carlo */
		nodes->clear();
		weights->clear();
		monte_carlo(nodes, weights, l, r);
		calculated_result = integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, Task_10::s0, Task_10::mu, Task_10::sigma, Task_10::T, Task_10::K);
		std::cout << calculated_result << std::endl;
		myfile<<calculate_relative_error(exact_result, calculated_result)<<"	";

		/* Trapezoidal rule */
		nodes->clear();
		weights->clear();
		trap_rule(nodes, weights, l);
		calculated_result = integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, Task_10::s0, Task_10::mu, Task_10::sigma, Task_10::T, Task_10::K);
		std::cout << calculated_result << std::endl;
		myfile<<calculate_relative_error(exact_result, calculated_result)<<"	";

		/* Clenshaw Curtis */
		nodes->clear();
		weights->clear();
		clenshaw_curtis(nodes, weights, l);
		calculated_result = integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, Task_10::s0, Task_10::mu, Task_10::sigma, Task_10::T, Task_10::K);
		std::cout << calculated_result << std::endl;
		myfile<<calculate_relative_error(exact_result, calculated_result)<<"	";

		/* Gauss Legendre */
		nodes->clear();
		weights->clear();
		gauss_legendre(nodes, weights, l);
		calculated_result = integrate_by_point_evaluation(call_option_integrand, nodes, weights, pow(2,l)-1, Task_10::s0, Task_10::mu, Task_10::sigma, Task_10::T, Task_10::K);
		std::cout << calculated_result << std::endl;
		myfile<<calculate_relative_error(exact_result, calculated_result)<<std::endl;
	}
    myfile.close();

	//frees memory
    gsl_rng_free(r);

	return 0;
}
