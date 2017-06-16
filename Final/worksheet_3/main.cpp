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
}

namespace Task_14{

	double T;
	int M;
}

namespace Task_15{
    bool write_in_file;
    bool use_trap_rule;
    bool use_bb;
	double simulation_result_rw;
	double simulation_result_bb;
	std::vector<std::vector<double> >* nodes;
	std::vector<std::vector<double> >* weights;
	int dimension_M;
	int l;
	double S0;
	double mu;
	double sigma;
	double T;
	double K;
}

namespace Task_16{
	double simulation_result_rw;
	double simulation_result_bb;
	int dimension_M;
	int l;
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

    Task_3::s0 = 10.;
    Task_3::r = 0.1;
    Task_3::T = 1.;
    Task_3::K = 10.;
    Task_3::sigma = 0.25;

	 /*Convergence plot for different N has to be inserted! */

//	std::cout << discrete_geometric_average_exact(Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma) << std::endl;

//	std::cout << discrete_geometric_average_simulation(rng, Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma, N) << std::endl;

//	std::cout << continuous_geometric_average_exact(Task_3::s0, Task_3::r, Task_3::T, Task_3::K, Task_3::sigma) << std::endl;
/*
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
*/
	/* Gauss-Legendre */
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

	/* Trapezoidal rule */
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

	/* Clenshaw Curtis */
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
	/*
	std::cout << "Van der Corput Sequence" << std::endl;
	std::vector<double> x = van_der_corput_sequence(3, 10, pow(10.,-12.));
	for(int i=0; i<10; i++)
		std::cout << x[i] << std::endl;

	std::cout << "Prime numbers" << std::endl;
	std::vector<int> prime_numbers = first_prime_numbers(20);
	for(int i=0; i<20; i++)
		std::cout << prime_numbers[i] << std::endl;
	*/

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

			std::cout << "QMC: " << integrate_by_point_evaluation_multivariate(function_task13, Task_13::N, d, Task_13::nodes, Task_13::weights, Task_13::gamma, d) << std::endl;


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
			std::cout << "N " << Task_13::N << std::endl;

			Task_13::nodes = new std::vector<std::vector<double>>((int)pow(8,d));
			Task_13::weights = new std::vector<double>((int)pow(8,d));
			full_grid_nodes_weights(Task_13::nodes, Task_13::weights, 7, d, trap_rule);

			Task_13::calculated_result = integrate_by_point_evaluation_multivariate(function_task13, (int)pow(7,d), Task_13::nodes, Task_13::weights, Task_13::gamma, d);
			std::cout << "Full grid trap rule: " << Task_13::calculated_result << std::endl;
			myfile << fabs(Task_13::exact_result - Task_13::calculated_result) << "	";

			// Full grid with Clenshaw Curtis
			std::cout << "N " << Task_13::N << std::endl;

			Task_13::nodes = new std::vector<std::vector<double>>((int)pow(8,d));
			Task_13::weights = new std::vector<double>((int)pow(8,d));
			full_grid_nodes_weights(Task_13::nodes, Task_13::weights, 7, d, clenshaw_curtis);

			Task_13::calculated_result = integrate_by_point_evaluation_multivariate(function_task13, (int)pow(7,d), Task_13::nodes, Task_13::weights, Task_13::gamma, d);
			std::cout << "Full grid Clenshaw Curtis: " << Task_13::calculated_result << std::endl;
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



  //////////////////////////////////////////////////////////
	//////////////////////////Task_14/////////////////////////
	//////////////////////////////////////////////////////////





	//////////////////////////////////////////////////////////
	//////////////////////////Task_15/////////////////////////
	//////////////////////////////////////////////////////////

	Task_15::S0 = 10;
	Task_15::T =1;
	Task_15::dimension_M = 16;
	Task_15::l = 1;
	Task_15::sigma = 0.25;
	Task_15::mu = 0.1;
	Task_15::K = 0;
	Task_15::nodes = new std::vector<std::vector<double> >(Task_15::dimension_M);
	Task_15::weights = new std::vector<std::vector<double> >(Task_15::dimension_M);
    Task_15::write_in_file = false;
    Task_15::use_trap_rule = false;
    Task_15::use_bb = true;

	double t = integrate_with_sparse_grid(asian_option_call_integrand,Task_15::dimension_M,Task_15::l,Task_15::nodes, Task_15::weights,Task_15::write_in_file,Task_15::use_trap_rule,Task_15::S0, Task_15::K,Task_15::sigma,Task_15::mu,Task_15::dimension_M,Task_15::T,Task_15::use_bb);
    	std::cout << "t= "<<t << '\n';



//////////////////////////////////////////////////////////
//////////////////////////Task_16/////////////////////////
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//////////////////////////Task_17/////////////////////////
//////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////
/*
	int d=10;
	int l=3;

	std::vector<std::vector<double> >* nodesv = new std::vector<std::vector<double> >(d);
	std::vector<std::vector<double> >* weightsv = new std::vector<std::vector<double> >(d);

	double final_value = integrate_with_sparse_grid(function_task13, d, l, nodesv, weightsv, true, false, 0.1, d);
	std::cout<<"final value"<<final_value<<std::endl;
	free(nodesv);
	free(weightsv);
	free(rng);
*/
}
