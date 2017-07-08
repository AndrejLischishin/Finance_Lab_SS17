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
#include "../header_files/multivariate_integration.hpp"
#include "../header_files/exotic_options.hpp"

namespace Task_1
{
	int plot_numbers;
	double s0;
	double r;
	double T;
	double K;
	double M;
	double sigma;
	double B;
	std::vector<double>* x;
	double payoff;
}

namespace Task_3
{
	double s0;
	double r;
	double T;
	double K;
	double sigma;
	double B;
	double fair_price;
}

namespace Task_4
{
	double s0;
	double K;
	double B;
	double T;
	double sigma;
	double r;
	std::vector<int> M(4);
	double exact_value;
	double calculated_result;
    std::vector<std::vector<double>>* nodes;
    std::vector<double>* weights_vec;
}

namespace Task_5
{
	int plot_numbers;
	double s0;
	double r;
	double T;
	double K;
	double M;
	double sigma;
	std::vector<double>* x;
	double payoff;
}

/**
 * Main function to run all exercises of worksheet 4.
 *
 * @param argc Integer argument for main
 * @param argv Char array argument for main
 *
 * @return Returns 0 if everything worked fine
 *
 */
int main(int argc, char* argv[])
{
	std::cout << "Prepared for worksheet 4!" << std::endl;

	gsl_rng* rng;

	//seeding
  	unsigned long seed = time(NULL);
  	//memory allocation
  	rng = gsl_rng_alloc(gsl_rng_mt19937);
  	gsl_rng_set(rng, seed);

	std::ofstream myfile;

	//////////////////////////////////////////////////////////
    //////////////////////////Task_1//////////////////////////
    //////////////////////////////////////////////////////////
	
	Task_1::plot_numbers = 100;
	Task_1::s0 = 10.0;
	Task_1::r = 0.02;
	Task_1::T = 1.0;
	Task_1::K = 10.0;
	Task_1::M = 2;
	Task_1::sigma = 0.2;
	Task_1::B = 8.0;

	Task_1::x = new std::vector<double>(2);

	myfile.open("output/plot_integrand_down_out_call.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	for(int i=1; i<Task_1::plot_numbers; i++)
	{
		for(int j=1; j<Task_1::plot_numbers; j++)
		{
			Task_1::x->clear();
			Task_1::x->push_back(i*Task_1::T/(Task_1::plot_numbers));
			Task_1::x->push_back(j*Task_1::T/(Task_1::plot_numbers));
			Task_1::payoff = payoff_discrete_down_out_call(Task_1::x, Task_1::s0, Task_1::r, Task_1::T, Task_1::M, Task_1::K, Task_1::sigma, Task_1::B);
			myfile << (*Task_1::x)[0] << "	" << (*Task_1::x)[1] << "	" << Task_1::payoff << std::endl;
		}
	}

	free(Task_1::x);

	myfile.close();

	//////////////////////////////////////////////////////////
    //////////////////////////Task_2//////////////////////////
    //////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////
    //////////////////////////Task_3//////////////////////////
    //////////////////////////////////////////////////////////

	Task_3::s0 = 10.0;
	Task_3::K = 10.0;
	Task_3::T = 1.0;
	Task_3::sigma = 0.2;
	Task_3::r = 0.02;

	myfile.open("output/fair_prices_down_out_call.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	for(int i=0; i<100; i++)
	{
		Task_3::B = 0.1*i;
		Task_3::fair_price = black_scholes_down_out_call(Task_3::s0, Task_3::K, Task_3::T, Task_3::sigma, Task_3::r, Task_3::B);
		myfile << Task_3::B << "	" << Task_3::fair_price << std::endl;
	}
	
	myfile.close();

	//////////////////////////////////////////////////////////
    //////////////////////////Task_4//////////////////////////
    //////////////////////////////////////////////////////////

	Task_4::s0 = 10.0;
	Task_4::K = 10.0;
	Task_4::B = 8.5;
	Task_4::T = 1.0;
	Task_4::sigma = 0.2;
	Task_4::r = 0.02;
	Task_4::M = {4,64,256,1024};

	std::cout << "Task 4" << std::endl;
	
	for(int n=10; n<=100000; n*=10)
	{
		Task_4::exact_value = black_scholes_down_out_call(Task_4::s0, Task_4::K, Task_4::T, Task_4::sigma, Task_4::r, Task_4::B);
		std::cout << "Exact result: " << Task_4::exact_value << std::endl;
		for(unsigned int m=0; m<Task_4::M.size(); m++)
		{
			Task_4::nodes = new std::vector<std::vector<double>>(n);
			if(Task_4::nodes==NULL)
				std::cout<<"Bad allocation"<<std::endl;
		
			Task_4::weights_vec = new std::vector<double>(n);
			monte_carlo_multivariate(Task_4::nodes, Task_4::weights_vec, n, Task_4::M[m], rng);

			Task_4::calculated_result = integrate_by_point_evaluation_multivariate(payoff_discrete_down_out_call, n, Task_4::nodes, Task_4::weights_vec, Task_4::s0, Task_4::r, Task_4::T, Task_4::M[m], Task_4::K, Task_4::sigma, Task_4::B);
			//myfile<<fabs(Task_4::exact_value-Task_4::calculated_result)<<" ";
			std::cout << Task_4::calculated_result << std::endl;
		}
	}
	
	free(Task_4::nodes);
	free(Task_4::weights_vec);

	//////////////////////////////////////////////////////////
    //////////////////////////Task_5//////////////////////////
    //////////////////////////////////////////////////////////

	Task_5::plot_numbers = 100;
	Task_5::s0 = 10.0;
	Task_5::r = 0.02;
	Task_5::T = 1.0;
	Task_5::K = 10.0;
	Task_5::M = 2;
	Task_5::sigma = 0.2;

	Task_5::x = new std::vector<double>(2);

	myfile.open("output/plot_integrand_lookback.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	for(int i=1; i<Task_5::plot_numbers; i++)
	{
		for(int j=1; j<Task_5::plot_numbers; j++)
		{
			Task_5::x->clear();
			Task_5::x->push_back(i*Task_5::T/(Task_5::plot_numbers));
			Task_5::x->push_back(j*Task_5::T/(Task_5::plot_numbers));
			Task_5::payoff = payoff_discrete_lookback(Task_5::x, Task_5::s0, Task_5::r, Task_5::T, Task_5::M, Task_5::K, Task_5::sigma);
			myfile << (*Task_5::x)[0] << "	" << (*Task_5::x)[1] << "	" << Task_5::payoff << std::endl;
		}
	}

	free(Task_5::x);

	myfile.close();

	//////////////////////////////////////////////////////////
    //////////////////////////Task_6//////////////////////////
    //////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////
    //////////////////////////Task_7//////////////////////////
    //////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////
    //////////////////////////Task_8//////////////////////////
    //////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////
    //////////////////////////Task_9//////////////////////////
    //////////////////////////////////////////////////////////

	return 0;
}
