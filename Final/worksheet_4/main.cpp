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

namespace Task_2 {
    
    int M;
    double s0;
    double K;
    double B;
    double T;
    double sigma;
    double r;
    
    double discretization;
    int N;
    
    
    std::vector<std::vector<double> >* nodes;
    std::vector<double>* weights_vec;
    std::vector<double> x;
    
    double reference_value;
    double calculated_result;
    bool use_bb;
    

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

namespace Task_6 {
    int M;
    double s0;
    double K;
    double T;
    double sigma;
    double r;
    
    double discretization;
    int N;

    
    std::vector<std::vector<double> >* nodes;
    std::vector<double>* weights_vec;
    std::vector<double> x;
    
    double reference_value;
    double calculated_result;
    bool use_bb;
}

namespace Task_7 {
    int N;
    int M;
    double s0;
    double K;
    double r;
    double sigma;
    double T;
    
    
    std::vector<std::vector<double> >* nodes;
    std::vector<double>* weights_vec;
    
    double geom_exact_value;
    double calculated_result;
    bool use_bb;
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
    
    Task_2::discretization = 64;
    Task_2::s0 = 10;
    Task_2::K = 10;
    Task_2::M = 64;
    Task_2::T = 1;
    Task_2::sigma = 0.2;
    Task_2::r = 0.02;
    Task_2::B = 8.5;
    
    Task_2::use_bb = true;
    
    myfile.open("output/task2.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

    for (int i = 0; i<Task_2::discretization; i++) {
        Task_2::x.push_back(i*Task_2::T/Task_2::discretization);
    }
    
    Task_2::reference_value = 0;
    //std::cout<<Task_2::reference_value<<std::endl;
    for (int i = 1; i<=4; i++) {
        
        
        Task_2::N = (int)pow(Task_2::M/2.,(double)i);
        myfile<<Task_2::N<<" ";
        
        //////////////////////////////////
        //////////////QMC_RW//////////////
        //////////////////////////////////
        
        Task_2::use_bb = false;

        
        Task_2::nodes = new std::vector<std::vector<double> >(Task_2::N);
        if (Task_2::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_2::weights_vec = new std::vector<double> (Task_2::N);
        if (Task_2::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        quasi_monte_carlo_multivariate(Task_2::nodes, Task_2::weights_vec, Task_2::N, Task_2::M);
        Task_2::calculated_result = integrate_by_point_evaluation_multivariate(barrier_integrand, Task_2::N, Task_2::nodes, Task_2::weights_vec, Task_2::s0,Task_2::K,Task_2::sigma,Task_2::r,Task_2::M,Task_2::T, Task_2::B, Task_2::use_bb);
        
        myfile<<Task_2::calculated_result<<" ";
        
        free(Task_2::nodes);
        free(Task_2::weights_vec);
        
        //////////////////////////////////
        //////////////QMC_BB//////////////
        //////////////////////////////////
        
        Task_2::use_bb = true;
        
        Task_2::nodes = new std::vector<std::vector<double> >(Task_2::N);
        if (Task_2::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_2::weights_vec = new std::vector<double> (Task_2::N);
        if (Task_2::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        quasi_monte_carlo_multivariate(Task_2::nodes, Task_2::weights_vec, Task_2::N, Task_2::M);
        Task_2::calculated_result = integrate_by_point_evaluation_multivariate(barrier_integrand, Task_2::N, Task_2::nodes, Task_2::weights_vec, Task_2::s0,Task_2::K,Task_2::sigma,Task_2::r,Task_2::M,Task_2::T, Task_2::B, Task_2::use_bb);
        
        myfile<<Task_2::calculated_result<<" ";
        
        free(Task_2::nodes);
        free(Task_2::weights_vec);
        
        //////////////////////////////////
        ///////////////MC_RW//////////////
        //////////////////////////////////
        
        Task_2::use_bb = false;
        
        Task_2::nodes = new std::vector<std::vector<double> >(Task_2::N);
        if (Task_2::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_2::weights_vec = new std::vector<double> (Task_2::N);
        if (Task_2::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        monte_carlo_multivariate(Task_2::nodes, Task_2::weights_vec, Task_2::N, Task_2::M,rng);
        Task_2::calculated_result = integrate_by_point_evaluation_multivariate(barrier_integrand, Task_2::N, Task_2::nodes, Task_2::weights_vec, Task_2::s0,Task_2::K,Task_2::sigma,Task_2::r,Task_2::M,Task_2::T,Task_2::B, Task_2::use_bb);
        
        myfile<<Task_2::calculated_result<<" ";
        
        free(Task_2::nodes);
        free(Task_2::weights_vec);
        
        //////////////////////////////////
        ///////////////MC_BB//////////////
        //////////////////////////////////
        
        Task_2::use_bb = true;
        
        Task_2::nodes = new std::vector<std::vector<double> >(Task_2::N);
        if (Task_2::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_2::weights_vec = new std::vector<double> (Task_2::N);
        if (Task_2::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        monte_carlo_multivariate(Task_2::nodes, Task_2::weights_vec, Task_2::N, Task_2::M,rng);
        Task_2::calculated_result = integrate_by_point_evaluation_multivariate(barrier_integrand, Task_2::N, Task_2::nodes, Task_2::weights_vec, Task_2::s0,Task_2::K,Task_2::sigma,Task_2::r,Task_2::M,Task_2::T,Task_2::B, Task_2::use_bb);
        
        myfile<<Task_2::calculated_result<<" ";
        myfile<<Task_2::reference_value<<std::endl;
        
        free(Task_2::nodes);
        free(Task_2::weights_vec);
        
    }
    myfile.close();

	//////////////////////////////////////////////////////////
    //////////////////////////Task_3//////////////////////////
    //////////////////////////////////////////////////////////


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

    Task_6::discretization = 64;
    Task_6::s0 = 10;
    Task_6::K = 10;
    Task_6::M = 64;
    Task_6::T = 1;
    Task_6::sigma = 0.2;
    Task_6::r = 0.02;
    
    Task_6::use_bb = true;
    
    myfile.open("output/task6.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    
    for (int i = 0; i<Task_6::discretization; i++) {
        Task_6::x.push_back(i*Task_6::T/Task_6::discretization);
    }
    
    Task_6::reference_value = payoff_discrete_lookback(Task_6::x, Task_6::s0, Task_6::r, Task_6::T, Task_6::M, Task_6::K, Task_6::sigma);
    //std::cout<<Task_6::reference_value<<std::endl;
    for (int i = 1; i<=4; i++) {
        
        
        Task_6::N = (int)pow(Task_6::M/2.,(double)i);
        myfile<<Task_6::N<<" ";
        
        //////////////////////////////////
        //////////////QMC/////////////////
        //////////////////////////////////
        
        Task_6::nodes = new std::vector<std::vector<double> >(Task_6::N);
        if (Task_6::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_6::weights_vec = new std::vector<double> (Task_6::N);
        if (Task_6::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        quasi_monte_carlo_multivariate(Task_6::nodes, Task_6::weights_vec, Task_6::N, Task_6::M);
        Task_6::calculated_result = integrate_by_point_evaluation_multivariate(lookback_call_integrand_fixed, Task_6::N, Task_6::nodes, Task_6::weights_vec, Task_6::s0,Task_6::K,Task_6::sigma,Task_6::r,Task_6::M,Task_6::T, Task_6::use_bb);
        
        myfile<<Task_6::calculated_result<<" ";
        
        free(Task_6::nodes);
        free(Task_6::weights_vec);
        
        //////////////////////////////////
        ///////////////MC/////////////////
        //////////////////////////////////
        
        Task_6::nodes = new std::vector<std::vector<double> >(Task_6::N);
        if (Task_6::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_6::weights_vec = new std::vector<double> (Task_6::N);
        if (Task_6::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        monte_carlo_multivariate(Task_6::nodes, Task_6::weights_vec, Task_6::N, Task_6::M,rng);
        Task_6::calculated_result = integrate_by_point_evaluation_multivariate(lookback_call_integrand_fixed, Task_6::N, Task_6::nodes, Task_6::weights_vec, Task_6::s0,Task_6::K,Task_6::sigma,Task_6::r,Task_6::M,Task_6::T, Task_6::use_bb);
        
        myfile<<Task_6::calculated_result<<" ";
        myfile<<Task_6::reference_value<<std::endl;
        
        free(Task_6::nodes);
        free(Task_6::weights_vec);
    
    }
    myfile.close();
	//////////////////////////////////////////////////////////
    //////////////////////////Task_7//////////////////////////
    //////////////////////////////////////////////////////////
    
    Task_7::s0 = 10;
    Task_7::T =1;
    Task_7::M = 64;
    Task_7::sigma = 0.25;
    Task_7::r = 0.1;
    Task_7::K = 10;
    
    
    Task_7::use_bb = true;
    
    myfile.open("task7.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    
    
    Task_7::geom_exact_value = discrete_geometric_average_exact(Task_7::s0, Task_7::r, Task_7::T, Task_7::M, Task_7::K, Task_7::sigma);
    
    for (int i = 1; i<=4; i++) {
        
        
        Task_7::N = (int)pow(Task_7::M/2.,(double)i);
        myfile<<Task_7::N<<" ";
        
        //////////////////////////////////
        //////////////QMC/////////////////
        //////////////////////////////////
        
        Task_7::nodes = new std::vector<std::vector<double> >(Task_7::N);
        if (Task_7::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_7::weights_vec = new std::vector<double> (Task_7::N);
        if (Task_7::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        quasi_monte_carlo_multivariate(Task_7::nodes, Task_7::weights_vec, Task_7::N, Task_7::M);
        Task_7::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand_arithmetic, Task_7::N, Task_7::nodes, Task_7::weights_vec, Task_7::s0,Task_7::K,Task_7::sigma,Task_7::r,Task_7::M,Task_7::T, Task_7::use_bb);
        
        myfile<<Task_7::calculated_result<<" ";
        
        free(Task_7::nodes);
        free(Task_7::weights_vec);
        
        //////////////////////////////////
        ///////QMC_control_variates///////
        //////////////////////////////////
        
        Task_7::nodes = new std::vector<std::vector<double> >(Task_7::N);
        if (Task_7::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_7::weights_vec = new std::vector<double> (Task_7::N);
        if (Task_7::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        quasi_monte_carlo_multivariate(Task_7::nodes, Task_7::weights_vec, Task_7::N, Task_7::M);
        Task_7::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand_arithmetic_control_variates, Task_7::N, Task_7::nodes, Task_7::weights_vec, Task_7::s0,Task_7::K,Task_7::sigma,Task_7::r, Task_7::M,Task_7::T, Task_7::use_bb);
        
        myfile<<Task_7::calculated_result+Task_7::geom_exact_value<<" ";
        
        free(Task_7::nodes);
        free(Task_7::weights_vec);
        
        //////////////////////////////////
        ///////////////MC/////////////////
        //////////////////////////////////
        
        Task_7::nodes = new std::vector<std::vector<double> >(Task_7::N);
        if (Task_7::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_7::weights_vec = new std::vector<double> (Task_7::N);
        if (Task_7::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        monte_carlo_multivariate(Task_7::nodes, Task_7::weights_vec, Task_7::N, Task_7::M, rng);
        Task_7::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand_arithmetic, Task_7::N, Task_7::nodes, Task_7::weights_vec, Task_7::s0,Task_7::K,Task_7::sigma,Task_7::r,Task_7::M,Task_7::T, Task_7::use_bb);
        
        myfile<<Task_7::calculated_result<<" ";
        
        free(Task_7::nodes);
        free(Task_7::weights_vec);
        
        //////////////////////////////////
        ////////MC_control_variates///////
        //////////////////////////////////
        
        Task_7::nodes = new std::vector<std::vector<double> >(Task_7::N);
        if (Task_7::nodes==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        Task_7::weights_vec = new std::vector<double> (Task_7::N);
        if (Task_7::weights_vec==NULL) {
            std::cout<<"Bad allocation task_7"<<std::endl;
        }
        
        monte_carlo_multivariate(Task_7::nodes, Task_7::weights_vec, Task_7::N, Task_7::M, rng);
        Task_7::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand_arithmetic_control_variates, Task_7::N, Task_7::nodes, Task_7::weights_vec, Task_7::s0,Task_7::K,Task_7::sigma,Task_7::r,Task_7::M,Task_7::T, Task_7::use_bb);
        
        myfile<<Task_7::calculated_result+Task_7::geom_exact_value<<std::endl;
        
        free(Task_7::nodes);
        free(Task_7::weights_vec);
        //std::cout<<"Lap"<<std::endl;
    }
    myfile.close();

	//////////////////////////////////////////////////////////
    //////////////////////////Task_8//////////////////////////
    //////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////
    //////////////////////////Task_9//////////////////////////
    //////////////////////////////////////////////////////////

    
    free(rng);
	return 0;
}