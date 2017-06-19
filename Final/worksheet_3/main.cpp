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

namespace Task_3
{
	double s0;
	double r;
	double T;
	int M;
	double K;
	double sigma;
}

namespace Task_4
{
	double s0;
	double r;
	double T;
	double K;
	double sigma;
}

namespace Task_5
{
	double s0;
	double r;
	double T;
	double K;
	double sigma;
	int M;
	int plot_numbers;
	std::vector<double> x;
	double payoff;
}

namespace Task_7
{
	int d;
	int n;
}

namespace Task_8
{
	int d;
	int l;
	std::vector<int> Nl;
	std::vector<double>* nodes;
	std::vector<double>* weights;
	std::vector<std::vector<double>> nodes_temp;
	std::vector<std::vector<double>> weights_temp;
	std::vector<int> ids;
}

namespace Task_9
{
	int d;
	int l;
	std::vector<int> Nl;
	std::vector<double>* nodes;
	std::vector<double>* weights;
	std::vector<std::vector<double>> nodes_temp;
	std::vector<int> ids;
}

namespace Task_12{
	int l;
}

namespace Task_13{
	int N;
	double gamma;
	std::vector<std::vector<double>>* nodes;
	std::vector<double>* weights;
	int max_l;
	double exact_result;
	double calculated_result;
    int N_full_grid;
}

namespace Task_14{

	double T;
	int M;
}

namespace Task_15{
  bool write_in_file;
  bool use_trap_rule;
  bool use_bb;
    double scale_factor;
	int max_level;
	double simulation_result_rw;
	double simulation_result_bb;
    double exact_value;
	std::vector<std::vector<double> >* nodes;
	std::vector<std::vector<double> >* weights;
	int dimension_M;
	long int N;

	double S0;
	double mu;
	double sigma;
	double T;
	double K;
}

namespace Task_16{
    bool write_in_file;
    bool use_trap_rule;
    bool use_bb;
    double scale_factor;
    double simulation_result_rw;
    double simulation_result_bb;
    double exact_value;
    double calculated_result;
    std::vector<std::vector<double> >* nodes;
    std::vector<std::vector<double> >* weights;
    std::vector<double>* weights_vec;
	int dimension_M;
    int max_level;
    long int N;

	double S0;
	double mu;
	double sigma;
	double T;
	double K;
}

namespace Task_17{
    bool write_in_file;
    bool use_trap_rule;
    bool use_bb;
    double scale_factor;
    double simulation_result_rw;
    double simulation_result_bb;
    double exact_value;
    double calculated_result;
    std::vector<std::vector<double> >* nodes;
    std::vector<std::vector<double> >* weights;
    std::vector<double>* weights_vec;
    int dimension_M;
    int max_level;
    int N;

    double S0;
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
 *
 */
int main(int argc, char* argv[])
{
	gsl_rng* rng;

	//seeding
  	unsigned long seed = time(NULL);
  	//memory allocation
  	rng = gsl_rng_alloc(gsl_rng_mt19937);
  	gsl_rng_set(rng, seed);

	std::ofstream myfile;

		//////////////////////////////////////////////////////////
    //////////////////////////Task_3//////////////////////////
    //////////////////////////////////////////////////////////
/*
    Task_3::s0 = 10.;
    Task_3::r = 0.1;
    Task_3::T = 1.;
    Task_3::K = 10.;
    Task_3::sigma = 0.25;

	 //Convergence plot for different N has to be inserted! 

//	std::cout << discrete_geometric_average_exact(Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma) << std::endl;

//	std::cout << discrete_geometric_average_simulation(rng, Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma, N) << std::endl;

//	std::cout << continuous_geometric_average_exact(Task_3::s0, Task_3::r, Task_3::T, Task_3::K, Task_3::sigma) << std::endl;

  Task_3::M = 10;
	double calculated_result;
	double exact_result;
	double absolute_error;

	myfile.open("output/error_task_3_M_10.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	for(int N=1; N<=1000000; N*=10)
	{
		calculated_result = discrete_geometric_average_simulation(rng, Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma, N);
		exact_result = discrete_geometric_average_exact(Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma);
		absolute_error = fabs(calculated_result-exact_result);
		myfile << N << " " << absolute_error << std::endl;
	}

	myfile.close();

	myfile.open("output/error_task_3_M_200.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	Task_3::M = 200;
	for(int N=1; N<=1000000; N*=10)
	{
		calculated_result = discrete_geometric_average_simulation(rng, Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma, N);
		exact_result = discrete_geometric_average_exact(Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma);
		absolute_error = fabs(calculated_result-exact_result);
		myfile << N << " " << absolute_error << std::endl;
	}

	myfile.close();

	//////////////////////////////////////////////////////////
	//////////////////////////Task_4//////////////////////////
	//////////////////////////////////////////////////////////

	Task_4::s0 = 10.;
	Task_4::r = 0.1;
	Task_4::T = 1.;
	Task_4::K = 10.;
	Task_4::sigma = 0.25;

	int M_max = pow(2,15);
	double continuous_result = continuous_geometric_average_exact(Task_4::s0, Task_4::r, Task_4::T, Task_4::K, Task_4::sigma);
	double discrete_result;
	double absolute_error;

	myfile.open("output/error_continuous_discrete.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	for(int M=1; M<=M_max; M*=2)
	{
		discrete_result = discrete_geometric_average_exact(Task_4::s0, Task_4::r, Task_4::T, M, Task_4::K, Task_4::sigma);
		absolute_error = fabs(continuous_result-discrete_result);
		myfile << M << " " << absolute_error << std::endl;
	}

	myfile.close();

	//////////////////////////////////////////////////////////
	//////////////////////////Task_5//////////////////////////
	//////////////////////////////////////////////////////////

	Task_5::plot_numbers = 100;
	Task_5::s0 = 10;
	Task_5::r = 0.1;
	Task_5::T = 1.0;
	Task_5::K = 10.0;
	Task_5::M = 2;
	Task_5::sigma = 0.25;

	myfile.open("output/plot_payoff_discrete_arithmetic_average.txt",std::ios::trunc);
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
			Task_5::payoff = payoff_discrete_arithmetic_average(Task_5::x, Task_5::s0, Task_5::r, Task_5::T, Task_5::M, Task_5::K, Task_5::sigma);
			myfile << Task_5::x[0] << "	" << Task_5::x[1] << "	" << Task_5::payoff << std::endl;
		}
	}

	myfile.close();


	//////////////////////////////////////////////////////////
	//////////////////////////Task_8//////////////////////////
	//////////////////////////////////////////////////////////
   
	Task_8::d=3;
	Task_8::l = 2;

  Task_8::nodes = new std::vector<double>;
  Task_8::weights = new std::vector<double>;
	gauss_legendre(Task_8::nodes, Task_8::weights, Task_8::l);

	for(int i=0; i<Task_8::d; i++)
    {
		Task_8::Nl.push_back((int)pow(2,Task_8::l)-1);
		std::vector<double> row;
		Task_8::nodes_temp.push_back(row);
		std::vector<double> row2;
		Task_8::weights_temp.push_back(row2);
		Task_8::ids.push_back(0);
    }

	for(int i=0; i<Task_8::d; i++)
	{
		for(int j=0; j<Task_8::Nl[i]; j++)
		{
			Task_8::nodes_temp[i].push_back((*Task_8::nodes)[j]);
			Task_8::weights_temp[i].push_back((*Task_8::weights)[j]);
		}
	}

	sum = 0.0;

	//tensor_product(0, Task_8::nodes_temp, Task_8::weights_temp, Task_8::d, Task_8::Nl, Task_8::ids, function_to_integrate);

	std::cout << sum << std::endl;

    



	//////////////////////////////////////////////////////////
	//////////////////////////Task_9//////////////////////////
	//////////////////////////////////////////////////////////

	Task_9::d = 2;
	Task_9::l = 5;

	// Gauss-Legendre
	Task_9::nodes = new std::vector<double>;
  Task_9::weights = new std::vector<double>;
	gauss_legendre(Task_9::nodes, Task_9::weights, Task_9::l);

	for(int i=0; i<Task_9::d; i++)
    {
		Task_9::Nl.push_back((int)pow(2,Task_9::l)-1);
		std::vector<double> row;
		Task_9::nodes_temp.push_back(row);
		Task_9::ids.push_back(0);
    }

	for(int i=0; i<Task_9::d; i++)
	{
		for(int j=0; j<Task_9::Nl[i]; j++)
		{
			Task_9::nodes_temp[i].push_back((*Task_9::nodes)[j]);
		}
	}

	myfile.open("output/quadrature_points_gauss_legendre.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
	//write_quadrature_points_to_file(myfile, 0, Task_9::nodes_temp, Task_9::d, Task_9::Nl, Task_9::ids);
	myfile.close();

	// Trapezoidal rule
	Task_9::nodes->clear();
    Task_9::weights->clear();
	Task_9::ids.clear();
	Task_9::Nl.clear();
	Task_9::nodes_temp.clear();
	trap_rule(Task_9::nodes, Task_9::weights, Task_9::l);

	for(int i=0; i<Task_9::d; i++)
    {
		Task_9::Nl.push_back((int)pow(2,Task_9::l)-1);
		std::vector<double> row;
		Task_9::nodes_temp.push_back(row);
		Task_9::ids.push_back(0);
    }

	for(int i=0; i<Task_9::d; i++)
	{
		for(int j=0; j<Task_9::Nl[i]; j++)
		{
			Task_9::nodes_temp[i].push_back((*Task_9::nodes)[j]);
		}
	}

	myfile.open("output/quadrature_points_trapezoidal_rule.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
	//write_quadrature_points_to_file(myfile, 0, Task_9::nodes_temp, Task_9::d, Task_9::Nl, Task_9::ids);
	myfile.close();

	// Clenshaw Curtis
	Task_9::nodes->clear();
    Task_9::weights->clear();
	Task_9::ids.clear();
	Task_9::Nl.clear();
	Task_9::nodes_temp.clear();
	clenshaw_curtis(Task_9::nodes, Task_9::weights, Task_9::l);

	for(int i=0; i<Task_9::d; i++)
    {
		Task_9::Nl.push_back((int)pow(2,Task_9::l)-1);
		std::vector<double> row;
		Task_9::nodes_temp.push_back(row);
		Task_9::ids.push_back(0);
    }

	for(int i=0; i<Task_9::d; i++)
	{
		for(int j=0; j<Task_9::Nl[i]; j++)
		{
			Task_9::nodes_temp[i].push_back((*Task_9::nodes)[j]);
		}
	}

	myfile.open("output/quadrature_points_clenshaw_curtis.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
	//write_quadrature_points_to_file(myfile, 0, Task_9::nodes_temp, Task_9::d, Task_9::Nl, Task_9::ids);
	myfile.close();

	//////////////////////////////////////////////////////////
	//////////////////////////Task_7//////////////////////////
	//////////////////////////////////////////////////////////
	
	std::cout << "Van der Corput Sequence" << std::endl;
	std::vector<double> x = van_der_corput_sequence(3, 10, pow(10.,-12.));
	for(int i=0; i<10; i++)
		std::cout << x[i] << std::endl;

	std::cout << "Prime numbers" << std::endl;
	std::vector<int> prime_numbers = first_prime_numbers(20);
	for(int i=0; i<20; i++)
		std::cout << prime_numbers[i] << std::endl;
	

	Task_7::d = 2;
	Task_7::n = 100;

	myfile.open("output/uniform_random_numbers.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	for(int i=0; i<Task_7::n; i++)
	{
		for(int j=0; j<Task_7::d; j++)
		{
			myfile << random_number_01() << "	";
		}
		myfile << std::endl;
	}

	myfile.close();

	myfile.open("output/halton_sequence.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

	//std::cout << "Halton sequence" << std::endl;
	std::vector<std::vector<double>> halton_sequence;
	halton_sequence = d_dimensional_halton_sequence(Task_7::d, Task_7::n);
	for(int i=0; i<Task_7::n; i++)
	{
		for(int j=0; j<Task_7::d; j++)
		{
			//std::cout << halton_sequence[i][j] << "  ";
			myfile << halton_sequence[i][j] << "	";
		}
		//std::cout << std::endl;
		myfile << std::endl;
	}

	myfile.close();

	//////////////////////////////////////////////////////////
	//////////////////////////Task_10/////////////////////////
	//////////////////////////////////////////////////////////

  	//////////////////////////////////////////////////////////
	//////////////////////////Task_11/////////////////////////
	//////////////////////////////////////////////////////////

  	//was done is already in pdf, latex

	//////////////////////////////////////////////////////////
	//////////////////////////Task_12/////////////////////////
	//////////////////////////////////////////////////////////

	myfile.open("output/number_of_points_SG_FG.txt",std::ios::trunc);
    if (!myfile.is_open()) {
       	std::cout<<"Error opening the file"<<std::endl;
    }

	Task_12::l = 4;

	for(int d=1; d<=10; d++)
	{
		int num_of_points_SG = (pow(2,Task_12::l)-1)*pow(Task_12::l,d-1);
		long int num_of_points_FG = pow((pow(2,Task_12::l)-1),d);

		std::cout << "num_of_points_SG: " << num_of_points_SG << std::endl;
		std::cout << "num_of_points_FG: " << num_of_points_FG << std::endl;

		myfile << d << "	" << num_of_points_SG << "	" << num_of_points_FG << std::endl;
	}

	myfile.close();

	*/
	
	//////////////////////////////////////////////////////////
	//////////////////////////Task_13/////////////////////////
	//////////////////////////////////////////////////////////


	Task_13::max_l = 4;

	Task_13::gamma = 0.1;

	for(int d=1; d<=8; d*=2)
	{
		Task_13::exact_result = function_task13_integral_exact_result(Task_13::gamma, d);
		std::cout << "Exact result: " << Task_13::exact_result << std::endl;

		myfile.open("output/testfunction_task13_error_d"+std::to_string(d)+".txt",std::ios::trunc);
    	if (!myfile.is_open()) {
        	std::cout<<"Error opening the file"<<std::endl;
    	}

		for(int l=1; l<=Task_13::max_l; l++)
		{
			Task_13::N = pow(2,l)*pow(l,d-1);

			myfile << Task_13::N << "	";

			// QMC
			Task_13::nodes = new std::vector<std::vector<double>>(Task_13::N);
			Task_13::weights = new std::vector<double>(Task_13::N);

			quasi_monte_carlo_multivariate(Task_13::nodes, Task_13::weights, Task_13::N, d);
			Task_13::calculated_result = integrate_by_point_evaluation_multivariate(function_task13, Task_13::N, Task_13::nodes, Task_13::weights, Task_13::gamma, d);
			std::cout << "QMC: " << Task_13::calculated_result << std::endl;
			myfile << fabs(Task_13::exact_result - Task_13::calculated_result) << "	";


			// MC
			Task_13::nodes = new std::vector<std::vector<double>>(Task_13::N);
			Task_13::weights = new std::vector<double>(Task_13::N);

			monte_carlo_multivariate(Task_13::nodes, Task_13::weights, Task_13::N, d, rng);
			Task_13::calculated_result = integrate_by_point_evaluation_multivariate(function_task13, Task_13::N, Task_13::nodes, Task_13::weights, Task_13::gamma, d);
			std::cout << "MC: " << Task_13::calculated_result << std::endl;
			myfile << fabs(Task_13::exact_result - Task_13::calculated_result) << "	";

			// Full grid with Trapezoidal rule
            Task_13::N_full_grid = (int)pow(Task_13::N, 1./d);

			Task_13::nodes = new std::vector<std::vector<double>>((int)pow(Task_13::N_full_grid,d+1));
			Task_13::weights = new std::vector<double>((int)pow(Task_13::N_full_grid,d+1));

			full_grid_nodes_weights(Task_13::nodes, Task_13::weights, Task_13::N_full_grid, d, trap_rule_absolute_number);
			Task_13::calculated_result = integrate_by_point_evaluation_multivariate(function_task13, (int)pow(Task_13::N_full_grid,d), Task_13::nodes, Task_13::weights, Task_13::gamma, d);
			std::cout << "Full grid trap rule: " << Task_13::calculated_result << std::endl;
			myfile << fabs(Task_13::exact_result - Task_13::calculated_result) << "	";

			// Full grid with Clenshaw Curtis
            Task_13::nodes = new std::vector<std::vector<double>>((int)pow(Task_13::N_full_grid,d+1));
            Task_13::weights = new std::vector<double>((int)pow(Task_13::N_full_grid,d+1));
            full_grid_nodes_weights(Task_13::nodes, Task_13::weights, Task_13::N_full_grid, d, clenshaw_curtis_absolute_number);
            
            Task_13::calculated_result = integrate_by_point_evaluation_multivariate(function_task13, (int)pow(Task_13::N_full_grid,d), Task_13::nodes, Task_13::weights, Task_13::gamma, d);
            std::cout << "Full grid clenshaw rule: " << Task_13::calculated_result << std::endl;
            myfile << fabs(Task_13::exact_result - Task_13::calculated_result) << "	";

			// Sparse Grid with Trapezoidal rule
			std::vector<std::vector<double>>* nodesv = new std::vector<std::vector<double>>(d);
			std::vector<std::vector<double>>* weightsv = new std::vector<std::vector<double>>(d);

			Task_13::calculated_result = integrate_with_sparse_grid(function_task13, d, l, nodesv, weightsv, false, true, Task_13::gamma, d);
			std::cout << "Sparse grid trap rule: " << Task_13::calculated_result << std::endl;
			myfile << fabs(Task_13::exact_result - Task_13::calculated_result) << "	";

			// Sparse Grid with Clenshaw Curtis
			nodesv = new std::vector<std::vector<double>>(d);
			weightsv = new std::vector<std::vector<double>>(d);

			Task_13::calculated_result = integrate_with_sparse_grid(function_task13, d, l, nodesv, weightsv, false, false, Task_13::gamma, d);
			std::cout << "Sparse grid Clenshaw Curtis: " << Task_13::calculated_result << std::endl;
			myfile << fabs(Task_13::exact_result - Task_13::calculated_result) << std::endl;
		}

		myfile.close();
	}

/*
	//////////////////////////////////////////////////////////
	//////////////////////////Task_14/////////////////////////
	//////////////////////////////////////////////////////////





	//////////////////////////////////////////////////////////
	//////////////////////////Task_15/////////////////////////
	//////////////////////////////////////////////////////////

	Task_15::S0 = 10;
	Task_15::T =1;
	Task_15::dimension_M = 16;

	Task_15::sigma = 0.25;
	Task_15::mu = 0.1;
	Task_15::K = 0;
    Task_15::scale_factor = exp(-Task_15::mu*Task_15::T);
	Task_15::nodes = new std::vector<std::vector<double> >(Task_15::dimension_M);
	Task_15::weights = new std::vector<std::vector<double> >(Task_15::dimension_M);
    Task_15::write_in_file = false;
    Task_15::use_trap_rule = false;
    Task_15::use_bb = true;
	Task_15::max_level = 4;
	Task_15::simulation_result_rw = 0;
	Task_15::simulation_result_bb = 0;

	myfile.open("output/testfunction_task15_error.txt",std::ios::trunc);
		if (!myfile.is_open()) {
				std::cout<<"Error opening the file"<<std::endl;
		}

    Task_15::exact_value = discrete_geometric_average_exact(Task_15::S0, Task_15::mu, Task_15::T, Task_15::dimension_M, Task_15::K, Task_15::sigma);
	for (int l = 1;l <= Task_15::max_level;l++) {
		Task_15::N = pow(2,l)*pow(l,Task_15::dimension_M-1);

		myfile << Task_15::N << "	";
        Task_15::simulation_result_rw = integrate_with_sparse_grid(asian_option_call_integrand,Task_15::dimension_M,l,Task_15::nodes, Task_15::weights,Task_15::write_in_file,Task_15::use_trap_rule,Task_15::S0, Task_15::K,Task_15::sigma,Task_15::mu,Task_15::dimension_M,Task_15::T,false);
		Task_15::simulation_result_bb = integrate_with_sparse_grid(asian_option_call_integrand,Task_15::dimension_M,l,Task_15::nodes, Task_15::weights,Task_15::write_in_file,Task_15::use_trap_rule,Task_15::S0, Task_15::K,Task_15::sigma,Task_15::mu,Task_15::dimension_M,Task_15::T,Task_15::use_bb);
        myfile << fabs(Task_15::scale_factor*Task_15::simulation_result_rw-Task_15::exact_value) << "	";
		myfile << fabs(Task_15::scale_factor*Task_15::simulation_result_bb-Task_15::exact_value) << std::endl;

	}
	myfile.close();


*/
//////////////////////////////////////////////////////////
//////////////////////////Task_16/////////////////////////
//////////////////////////////////////////////////////////
    
    Task_16::S0 = 10;
    Task_16::T =1;
    Task_16::dimension_M = 8;
    
    Task_16::sigma = 0.25;
    Task_16::mu = 0.1;
    Task_16::K = 0;
    Task_16::scale_factor = exp(-Task_16::mu*Task_16::T);

    Task_16::write_in_file = false;
    Task_16::use_trap_rule = false;
    Task_16::use_bb = true;
    Task_16::max_level = 4;
    Task_16::simulation_result_rw = 0;
    Task_16::simulation_result_bb = 0;
    
    myfile.open("output/testfunction_task16_error.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    
    
    Task_16::exact_value = discrete_geometric_average_exact(Task_16::S0, Task_16::mu, Task_16::T, Task_16::dimension_M, Task_16::K, Task_16::sigma);
    
    for (int l = 1; l<=Task_16::max_level; l++) {
        
        Task_16::N = pow(2,l)*pow(l,Task_16::dimension_M-1);
        myfile<<Task_16::N<<" ";
        
        //////////////////////////////////
        //////////////QMC/////////////////
        //////////////////////////////////
        Task_16::nodes = new std::vector<std::vector<double> >(Task_16::N);
        Task_16::weights_vec = new std::vector<double> (Task_16::N);
        //QMC without Brownian Bridge
        quasi_monte_carlo_multivariate(Task_16::nodes, Task_16::weights_vec, (int)Task_16::N, Task_16::dimension_M);
        Task_16::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)Task_16::N, Task_16::nodes, Task_16::weights_vec, Task_16::S0,Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,false);
        std::cout<<Task_16::calculated_result<<std::endl;
        myfile<<fabs(Task_16::exact_value-Task_16::scale_factor*Task_16::calculated_result)<<" ";
        //QMC with Brownian Bridge
        Task_16::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)Task_16::N, Task_16::nodes, Task_16::weights_vec, Task_16::S0,Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,true);
        std::cout<<Task_16::calculated_result<<std::endl;
        myfile<<fabs(Task_16::exact_value-Task_16::scale_factor*Task_16::calculated_result)<<" ";
        free(Task_16::nodes);
        free(Task_16::weights_vec);
        //////////////////////////////////
        //////////////MC//////////////////
        //////////////////////////////////
        Task_16::nodes = new std::vector<std::vector<double> >(Task_16::N);
        Task_16::weights_vec = new std::vector<double> (Task_16::N);
        //MC without Brownian Bridge
        monte_carlo_multivariate(Task_16::nodes, Task_16::weights_vec, (int)Task_16::N, Task_16::dimension_M,rng);
        Task_16::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)Task_16::N, Task_16::nodes, Task_16::weights_vec, Task_16::S0,Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,false);
        std::cout<<Task_16::calculated_result<<std::endl;
        myfile<<fabs(Task_16::exact_value-Task_16::scale_factor*Task_16::calculated_result)<<" ";
        //MC with Brownian Bridge
        Task_16::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)Task_16::N, Task_16::nodes, Task_16::weights_vec, Task_16::S0,Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,true);
        myfile<<fabs(Task_16::exact_value-Task_16::scale_factor*Task_16::calculated_result)<<" ";
        std::cout<<Task_16::calculated_result<<std::endl;
        free(Task_16::nodes);
        free(Task_16::weights_vec);
        //////////////////////////////////
        /////Full_grid_with_TRAPEZOIDAl///
        //////////////////////////////////
        //without_Brownian_Bridge
        Task_16::nodes = new std::vector<std::vector<double> >((int)pow((pow(2,l)-1),Task_16::dimension_M+1));
        Task_16::weights_vec = new std::vector<double> ((int)pow((pow(2,l)-1),Task_16::dimension_M+1));
        
        full_grid_nodes_weights(Task_16::nodes, Task_16::weights_vec, (pow(2,l)-1), Task_16::dimension_M, trap_rule_absolute_number);
        
        Task_16::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)pow((pow(2,l)-1),Task_16::dimension_M), Task_16::nodes, Task_16::weights_vec, Task_16::S0,Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,false);
        std::cout<<" asdf "<<Task_16::calculated_result<<std::endl;
        myfile<<fabs(Task_16::exact_value-Task_16::scale_factor*Task_16::calculated_result)<<" ";
        //with_Brownian_Bridge
        Task_16::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)pow((pow(2,l)-1),Task_16::dimension_M), Task_16::nodes, Task_16::weights_vec, Task_16::S0,Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,true);
        std::cout<<Task_16::calculated_result<<std::endl;
        myfile<<fabs(Task_16::exact_value-Task_16::scale_factor*Task_16::calculated_result)<<" ";
        free(Task_16::nodes);
        free(Task_16::weights_vec);
        //////////////////////////////////
        /////Full_grid_with_CC////////////
        //////////////////////////////////
        //without_Brownian_Bridge
        Task_16::nodes = new std::vector<std::vector<double> >((int)pow((pow(2,l)),Task_16::dimension_M));
        Task_16::weights_vec = new std::vector<double> ((int)pow((pow(2,l)),Task_16::dimension_M));
        full_grid_nodes_weights(Task_16::nodes, Task_16::weights_vec, (pow(2,l)-1), Task_16::dimension_M, clenshaw_curtis_absolute_number);
        
        Task_16::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)pow((pow(2,l)-1),Task_16::dimension_M), Task_16::nodes, Task_16::weights_vec, Task_16::S0,Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,false);
        std::cout<<Task_16::calculated_result<<std::endl;
        myfile<<fabs(Task_16::exact_value-Task_16::scale_factor*Task_16::calculated_result)<<" ";
        //with_Brownian_Bridge
        Task_16::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)pow((pow(2,l)-1),Task_16::dimension_M), Task_16::nodes, Task_16::weights_vec, Task_16::S0,Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,true);
        std::cout<<Task_16::calculated_result<<std::endl;
        myfile<<fabs(Task_16::exact_value-Task_16::scale_factor*Task_16::calculated_result)<<" ";
        
        free(Task_16::nodes);
        free(Task_16::weights_vec);
        //////////////////////////////////
        ///Sparse_grid_with_TRAPEZOIDAl///
        //////////////////////////////////
        
        Task_16::nodes = new std::vector<std::vector<double> >(Task_16::dimension_M);
        Task_16::weights = new std::vector<std::vector<double> > (Task_16::dimension_M);
        
        //without_Brownian_Bridge
        Task_16::simulation_result_rw = integrate_with_sparse_grid(asian_option_call_integrand,Task_16::dimension_M,l,Task_16::nodes, Task_16::weights,Task_16::write_in_file,true,Task_16::S0, Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,false);
        std::cout<<"sG "<<Task_16::simulation_result_rw<<std::endl;
        myfile << fabs(Task_16::scale_factor*Task_16::simulation_result_rw-Task_16::exact_value)<<" ";
        //with_Brownian_Bridge
        Task_16::simulation_result_bb = integrate_with_sparse_grid(asian_option_call_integrand,Task_16::dimension_M,l,Task_16::nodes, Task_16::weights,Task_16::write_in_file,true,Task_16::S0, Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,true);
        std::cout<<Task_16::simulation_result_bb<<std::endl;
        myfile << fabs(Task_16::scale_factor*Task_16::simulation_result_bb-Task_16::exact_value)<<" ";
        
        free(Task_16::nodes);
        free(Task_16::weights);
        //////////////////////////////////
        ///Sparse_grid_with_CC////////////
        //////////////////////////////////
        
        Task_16::nodes = new std::vector<std::vector<double> >(Task_16::dimension_M);
        Task_16::weights = new std::vector<std::vector<double> > (Task_16::dimension_M);
        
        //without_Brownian_Bridge
        Task_16::simulation_result_rw = integrate_with_sparse_grid(asian_option_call_integrand,Task_16::dimension_M,l,Task_16::nodes, Task_16::weights,Task_16::write_in_file,false,Task_16::S0, Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,false);
        std::cout<<Task_16::simulation_result_rw<<std::endl;
        myfile << fabs(Task_16::scale_factor*Task_16::simulation_result_rw-Task_16::exact_value) <<" ";
        //with_Brownian_Bridge
        Task_16::simulation_result_bb = integrate_with_sparse_grid(asian_option_call_integrand,Task_16::dimension_M,l,Task_16::nodes, Task_16::weights,Task_16::write_in_file,false,Task_16::S0, Task_16::K,Task_16::sigma,Task_16::mu,Task_16::dimension_M,Task_16::T,true);
        std::cout<<Task_16::simulation_result_bb<<std::endl;
        myfile << fabs(Task_16::scale_factor*Task_16::simulation_result_bb-Task_16::exact_value)<<std::endl;
        
        free(Task_16::nodes);
        free(Task_16::weights);
        
    }
    

    myfile.close();

//////////////////////////////////////////////////////////
//////////////////////////Task_17/////////////////////////
//////////////////////////////////////////////////////////
/*

    Task_17::S0 = 10;
    Task_17::T =1;
    Task_17::dimension_M = 64;

    Task_17::sigma = 0.25;
    Task_17::mu = 0.1;
    Task_17::K = 10;
    Task_17::scale_factor = exp(-Task_17::mu*Task_17::T);
    Task_17::write_in_file = false;
    Task_17::use_trap_rule = false;
    Task_17::use_bb = true;
    Task_17::max_level = 4;
    Task_17::simulation_result_rw = 0;
    Task_17::simulation_result_bb = 0;

    myfile.open("output/testfunction_task17_error.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }


    Task_17::exact_value = discrete_geometric_average_exact(Task_17::S0, Task_17::mu, Task_17::T, Task_17::dimension_M, Task_17::K, Task_17::sigma);

    for (int l = 1; l<=Task_17::max_level; l++) {

        Task_17::N = Task_17::dimension_M*l*(int)pow(10, (double)l);
        myfile<<Task_17::N<<" ";

        //////////////////////////////////
        //////////////QMC/////////////////
        //////////////////////////////////


        Task_17::nodes = new std::vector<std::vector<double> >(Task_17::N);
        if (Task_17::nodes==NULL) {
            std::cout<<"Bad allocation"<<std::endl;
        }
        Task_17::weights_vec = new std::vector<double> (Task_17::N);
        //QMC without Brownian Bridge
        quasi_monte_carlo_multivariate(Task_17::nodes, Task_17::weights_vec, (int)Task_17::N, Task_17::dimension_M);

        Task_17::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)Task_17::N, Task_17::nodes, Task_17::weights_vec, Task_17::S0,Task_17::K,Task_17::sigma,Task_17::mu,Task_17::dimension_M,Task_17::T,false);
        myfile<<fabs(Task_17::exact_value-Task_17::scale_factor*Task_17::calculated_result)<<" ";
        //QMC with Brownian Bridge
        Task_17::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)Task_17::N, Task_17::nodes, Task_17::weights_vec, Task_17::S0,Task_17::K,Task_17::sigma,Task_17::mu,Task_17::dimension_M,Task_17::T,true);
        myfile<<fabs(Task_17::exact_value-Task_17::scale_factor*Task_17::calculated_result)<<" ";
        free(Task_17::nodes);
        free(Task_17::weights_vec);

        //////////////////////////////////
        //////////////MC//////////////////
        //////////////////////////////////
        Task_17::nodes = new std::vector<std::vector<double> >(Task_17::N);
        Task_17::weights_vec = new std::vector<double> (Task_17::N);
        //MC without Brownian Bridge
        monte_carlo_multivariate(Task_17::nodes, Task_17::weights_vec, (int)Task_17::N, Task_17::dimension_M,rng);
        Task_17::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)Task_17::N, Task_17::nodes, Task_17::weights_vec, Task_17::S0,Task_17::K,Task_17::sigma,Task_17::mu,Task_17::dimension_M,Task_17::T,false);
        myfile<<fabs(Task_17::exact_value-Task_17::scale_factor*Task_17::calculated_result)<<" ";
        //MC with Brownian Bridge
        Task_17::calculated_result = integrate_by_point_evaluation_multivariate(asian_option_call_integrand, (int)Task_17::N, Task_17::nodes, Task_17::weights_vec, Task_17::S0,Task_17::K,Task_17::sigma,Task_17::mu,Task_17::dimension_M,Task_17::T,true);
        myfile<<fabs(Task_17::exact_value-Task_17::scale_factor*Task_17::calculated_result)<<std::endl;
        free(Task_17::nodes);
        free(Task_17::weights_vec);
    }


    myfile.close();

 */
    free(rng);

/////////////////////////////////////////////////////////
}
