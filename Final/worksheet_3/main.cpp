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


double sum;

/**
 * Put the functions into header files later on!
 */
//////////////////////////////////////////////////////////
//////////////////////////Task_2//////////////////////////
//////////////////////////////////////////////////////////

double discrete_geometric_average_exact(double s0, double r, double T, int M, double K, double sigma)
{
	double delta_t = T/(double)M;
	double T1 = T-((M*(M-1.)*(4.*M+1.))/(6.*M*M))*delta_t;
	double T2 = T-((M-1.)/2.)*delta_t;
	double A = exp(-r*(T-T2)-(sigma*sigma*(T2-T1))/2);
	double d = (log(s0/K)+(r-0.5*sigma*sigma)*T2)/(sigma*sqrt(T1));

	return s0*A*normal_cdf(d+sigma*sqrt(T1))-K*exp(-r*T)*normal_cdf(d);
}

double continuous_geometric_average_exact(double s0, double r, double T, double K, double sigma)
{
	double d = (log(s0/K)+0.5*(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T/3.));

	return s0*exp(-0.5*(r+(sigma*sigma)/6.)*T)*normal_cdf(d+sigma*sqrt(T/3.))-K*exp(-r*T)*normal_cdf(d);
}

double discrete_geometric_average_simulation(gsl_rng* rng, double s0, double r, double T, int M, double K, double sigma, int N)
{
	double delta_t = T/(double)M;
	double result = 0.0;
	std::vector<double>* w;
	std::vector<double>* s;

	for(int i=0; i<N; i++)
	{
		/* Simulation of brownian motion */
		w = wiener_process(rng, T, delta_t);
		s = brownian_motion(rng, T, delta_t, w, s0, r, sigma);

		double product = 1.0;
		for(int j=0; j<M; j++)
		{
			product *= (*s)[j];
		}
		product = pow(product, 1./(double)M)-K;
		if(product > 0)
			result += product;
	}
	result = result/(double)N;
	return result;
}


//////////////////////////////////////////////////////////
//////////////////////////Task_8//////////////////////////
//////////////////////////////////////////////////////////
/* Nl has all the different N_l_1,....,N_l_d */
double function_to_integrate(std::vector<double> x)
{
	return x[0]+x[1]+x[2];
}


template<typename... Args>
void tensor_product(int iteration, std::vector<std::vector<double>> nodes_temp, std::vector<std::vector<double>> weights_temp, int d, std::vector<int> Nl, std::vector<int> ids, double (*function_to_integrate)(std::vector<double> x, Args... rest), Args... rest)
{
	if(iteration==d)
	{ 
		std::vector<double> x;
		x.clear();
		double prod = 1.0;
		//std::cout << "Test" << std::endl;
		for(int i=0; i<d; i++)
		{
			x.push_back(nodes_temp[i][ids[i]]);
			prod *= weights_temp[i][ids[i]];
			//std::cout << nodes_temp[i][ids[i]] << " ";
		}
		sum += prod*function_to_integrate(x, rest...);
		//std::cout << sum << std::endl;
	}
	else
	{
		for(int k=0; k<Nl[iteration]; k++)
		{
			ids[iteration] = k;
			tensor_product(iteration+1, nodes_temp, weights_temp, d, Nl, ids, function_to_integrate, rest...);		
		}
	}
}

//////////////////////////////////////////////////////////
//////////////////////////Task_9//////////////////////////
//////////////////////////////////////////////////////////
template<typename... Args>
void write_quadrature_points_to_file(std::ofstream& myfile, int iteration, std::vector<std::vector<double>> nodes_temp, int d, std::vector<int> Nl, std::vector<int> ids)
{
	if(iteration==d)
	{
		for(int i=0; i<d; i++)
		{
			myfile << nodes_temp[i][ids[i]] << "	";
			//std::cout << nodes_temp[i][ids[i]] << "	";
		}
		myfile << std::endl;
		//std::cout << std::endl;
	}
	else
	{
		for(int k=0; k<Nl[iteration]; k++)
		{
			ids[iteration] = k;
			write_quadrature_points_to_file(myfile, iteration+1, nodes_temp, d, Nl, ids);		
		}
	}
}





namespace Task_3
{
	double s0;
	double r;
	double T;
	int M;
	double K;
	double sigma;
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

	std::cout << "Prepared everything for worksheet 3." << std::endl;

	//////////////////////////////////////////////////////////
    //////////////////////////Task_3//////////////////////////
    //////////////////////////////////////////////////////////

	Task_3::s0 = 10.;
    Task_3::r = 0.1;
    Task_3::T = 1.;										
    Task_3::M = 10;
    Task_3::K = 10.;	
	Task_3::sigma = 0.25;
	
	int N = 1000000;

	/* Convergence plot for different N has to be inserted! */

	std::cout << discrete_geometric_average_exact(Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma) << std::endl;
	
	std::cout << discrete_geometric_average_simulation(rng, Task_3::s0, Task_3::r, Task_3::T, Task_3::M, Task_3::K, Task_3::sigma, N) << std::endl;

	std::cout << continuous_geometric_average_exact(Task_3::s0, Task_3::r, Task_3::T, Task_3::K, Task_3::sigma) << std::endl;


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

	tensor_product(0, Task_8::nodes_temp, Task_8::weights_temp, Task_8::d, Task_8::Nl, Task_8::ids, function_to_integrate);	

	std::cout << sum << std::endl;

	//////////////////////////////////////////////////////////
	//////////////////////////Task_9//////////////////////////
	//////////////////////////////////////////////////////////

	Task_9::d = 2;
	Task_9::l = 5;

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
	
	std::ofstream myfile;
	myfile.open("output/quadrature_points_gauss_legendre.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
	write_quadrature_points_to_file(myfile, 0, Task_9::nodes_temp, Task_9::d, Task_9::Nl, Task_9::ids);
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
	write_quadrature_points_to_file(myfile, 0, Task_9::nodes_temp, Task_9::d, Task_9::Nl, Task_9::ids);
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
	write_quadrature_points_to_file(myfile, 0, Task_9::nodes_temp, Task_9::d, Task_9::Nl, Task_9::ids);
	myfile.close();

	return 0;
}
