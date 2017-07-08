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
	std::vector<double> x;
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

namespace Task_5
{
	int plot_numbers;
	double s0;
	double r;
	double T;
	double K;
	double M;
	double sigma;
	std::vector<double> x;
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

	myfile.open("output/plot_integrand_down_out_call.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	for(int i=1; i<Task_1::plot_numbers; i++)
	{
		for(int j=1; j<Task_1::plot_numbers; j++)
		{
			Task_1::x.clear();
			Task_1::x.push_back(i*Task_1::T/(Task_1::plot_numbers));
			Task_1::x.push_back(j*Task_1::T/(Task_1::plot_numbers));
			Task_1::payoff = payoff_discrete_down_out_call(Task_1::x, Task_1::s0, Task_1::r, Task_1::T, Task_1::M, Task_1::K, Task_1::sigma, Task_1::B);
			myfile << Task_1::x[0] << "	" << Task_1::x[1] << "	" << Task_1::payoff << std::endl;
		}
	}

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

	myfile.open("output/plot_integrand_lookback.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	for(int i=1; i<Task_5::plot_numbers; i++)
	{
		for(int j=1; j<Task_5::plot_numbers; j++)
		{
			Task_5::x.clear();
			Task_5::x.push_back(i*Task_5::T/(Task_5::plot_numbers));
			Task_5::x.push_back(j*Task_5::T/(Task_5::plot_numbers));
			Task_5::payoff = payoff_discrete_lookback(Task_5::x, Task_5::s0, Task_5::r, Task_5::T, Task_5::M, Task_5::K, Task_5::sigma);
			myfile << Task_5::x[0] << "	" << Task_5::x[1] << "	" << Task_5::payoff << std::endl;
		}
	}

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
