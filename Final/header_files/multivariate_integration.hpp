
#ifndef multivariate_integration_hpp
#define multivariate_integration_hpp

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include "integration_functions.hpp"
#include "random_functions.hpp"
#include "simulation_functions.hpp"
//trap_rule
void trap_rule_weights(std::vector<double>* weights, int l);
void trap_rule_nodes(std::vector<double> *nodes, int l);
//clenshaw curtis
void clenshaw_curtis_weights(std::vector<double>* weights, int l);
void clenshaw_curtis_nodes(std::vector<double> *nodes, int l);
//sparse grid
void sparse_grid_nodes(int d, int product, std::vector<int>* allvec, std::vector<std::vector<double> >* nodesv);
//utility
int enumeration(std::vector<int>* k, int d, int l, std::vector<int>* diag);
void loop(std::vector<int>* vec, std::vector<int>* klevel, int d, std::vector<int>* finalvec);
//finance functions
double discrete_geometric_average_exact(double s0, double r, double T, int M, double K, double sigma);
double continuous_geometric_average_exact(double s0, double r, double T, double K, double sigma);
double discrete_geometric_average_simulation(gsl_rng* rng, double s0, double r, double T, int M, double K, double sigma, int N);
double asian_option_call_integrand(std::vector<double> x,double S0,double K, double sigma, double mu,int M, double T, bool bb_not_rw);
//sequences
std::vector<double> van_der_corput_sequence(int p, int n, double epsilon);
bool is_prime(int number);
std::vector<int> first_prime_numbers(int n);
std::vector<std::vector<double>> d_dimensional_halton_sequence(int d, int n);
std::vector<double> van_der_corput_sequence(int p, int n, double epsilon);
//brownian bridge
std::vector<double>* brownian_bridge_level(std::vector<double>* prev_level, double T, double level, gsl_rng* r);
std::vector<double>* brownian_bridge(gsl_rng* r, double T, int M);
//task 13
double function_task13(std::vector<double>* x, double gamma, int d);
double function_task13_sparse(std::vector<double>* x, double gamma, int d);
//utility
void write_quadrature_points_to_file(std::ofstream& myfile, int iteration, std::vector<std::vector<double>> nodes_temp, int d, std::vector<int> Nl, std::vector<int> ids);
double function_to_integrate(std::vector<double> x);
//Monte Carlo
void monte_carlo_multivariate(std::vector<std::vector<double>>* nodes, std::vector<double>* weights, int l, int d, gsl_rng* r);
void quasi_monte_carlo_multivariate(std::vector<std::vector<double>>* nodes, std::vector<double>* weights, int l, int d);

//////////////////////////////////////////////////
///////////tensor_product/////////////////////////
//////////////////////////////////////////////////
//double sum;

/*template<typename... Args>
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
 */
//////////////////////////////////////////////////
///integrate_with_sparse_grid/////////////////////
//////////////////////////////////////////////////
template<typename... Args>
double integrate_with_sparse_grid(double (*multifunction_to_integrate)(std::vector<double>* x, Args... rest),int d, int l, std::vector<std::vector<double>>* nodes,std::vector<std::vector<double>>* weights, bool write_in_file, bool use_trap_rule, Args... rest){
	if(write_in_file==true) {
		FILE *fp;
	fp = fopen("stuetzstellen", "w");
		fprintf(fp,"");

	fclose(fp);

	};
	int maxlevel = (int)pow(2,l)-1;
//	int* klevel = new int[d];
	std::vector<int>* klevel = new std::vector<int>(d);
//	int* k = new int[d];
	std::vector<int>* k = new std::vector<int>(d);
//	int* k1 = new int[d];
	std::vector<int>* k1 = new std::vector<int>(d);
//	int* vec = new int[d];
	std::vector<int>* vec = new std::vector<int>(d);
	int sum = 1;
	std::vector<int>* diag = new std::vector<int>;

	int K;
	int J;
	double product_w=1;
	int product = 1;
	double final_value=0;

	for(int i=0; i<d; i++){
		(*vec)[i]=1;
		(*k1)[i]=1;
		sum = sum * maxlevel;
	}

	std::vector<int>* allvec = new std::vector<int>;

	int sz = enumeration(k1,d,l,diag);

	for(K=0; K<sz; K++)
	{
	// von diag auf k uebertragen
		for(int i=0; i<d; i++){

			(*k)[i]=(*diag)[i+K*d]+1;
		}

	// klevel
		for(int i=0; i<d; i++){
			(*klevel)[i]=pow(2,(*k)[i])-1;
		}

	// vec auf 1 setzen
		for(int i=0; i<d; i++){
			(*vec)[i]=1;
		}


		// Tensorprodukt

		loop(vec, klevel, d, allvec);

		// produkt berechnen
		product = 1;
		for(int i=0; i<d; i++){
			product = product * (*klevel)[i];
		}

		if(use_trap_rule== true){
			for(int i=0; i<d; i++){
				trap_rule_weights(&(*weights)[i], (*k)[i]);
				trap_rule_nodes(&(*nodes)[i], (*k)[i]);
			}
		}
		else{
			for(int i=0; i<d; i++){
				clenshaw_curtis_weights(&(*weights)[i], (*k)[i]);
				clenshaw_curtis_nodes(&(*nodes)[i], (*k)[i]);
			}
		}



		for(int i=0; i<product*d; i=i+d){
			for(int j=0; j<d; j++){
				product_w = product_w*(*weights)[j][(*allvec)[i+j]];

				J =j;
			}
			std::vector<double>* point = new std::vector<double>;
			point->clear();
			for(int h=0; h<d; h++){
				point->push_back((*nodes)[h][(*allvec)[i+h]]);
			}
			final_value = final_value + product_w*(multifunction_to_integrate(point, rest...));
			product_w = 1;


	if(write_in_file==true) sparse_grid_nodes(d,product,allvec,nodes);

		}


	}
	free(allvec);
	free(k1);
	free(k);
	free(vec);
	free(diag);

	return final_value;

}
//////////////////////////////////////////////////
///integrate_by_point_evaluation_multivariate/////
//////////////////////////////////////////////////

template<typename... Args>
 double integrate_by_point_evaluation_multivariate(double (*function)(std::vector<double>* x, Args... rest),int n, int d, std::vector<std::vector<double>>* nodes, std::vector<double>* weights, Args... rest)
 {
   double result = 0.0;
   for(int i=0; i<pow(2,n)-1; i++)
   {
     result += (*weights)[i]*function(&(*nodes)[i], rest...);
   }

   return result;
 }


#endif /* multivariate_integration_hpp */
