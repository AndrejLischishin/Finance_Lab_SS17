
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

double testfunction(std::vector<double> x);
void trap_rule_s(std::vector<double>* weights, int l);
void trap_rulen(std::vector<double> *nodes, int l);
int enumeration(int *k, int d, int l, std::vector<int>* diag);
void loop(int* vec, int* klevel, int d, int* finalvec);
void sparse_grid_nodes(int d, int product, int* allvec, double** nodes, std::vector<std::vector<double> > nodesv);
void sparse_grid_weights(int d, int product, int* allvec, double** weights, std::vector<std::vector<double> > weightsv);
double discrete_geometric_average_exact(double s0, double r, double T, int M, double K, double sigma);
double continuous_geometric_average_exact(double s0, double r, double T, double K, double sigma);
double discrete_geometric_average_simulation(gsl_rng* rng, double s0, double r, double T, int M, double K, double sigma, int N);
double asian_option_call_integrand(std::vector<double> x,double S0,double K, double sigma, double mu,int M, double T, bool bb_not_rw);
double function_to_integrate(std::vector<double> x);
std::vector<double> van_der_corput_sequence(int p, int n, double epsilon);
bool is_prime(int number);
std::vector<int> first_prime_numbers(int n);
std::vector<std::vector<double>> d_dimensional_halton_sequence(int d, int n);
std::vector<double>* brownian_bridge_level(std::vector<double>* prev_level, double T, double level, gsl_rng* r);
std::vector<double>* brownian_bridge(gsl_rng* r, double T, int M);
double function_task13(std::vector<double>* x, double gamma, int d);
double function_task13_sparse(std::vector<double> x, double gamma, int d);

//template<typename... Args>
//void tensor_product(int iteration, std::vector<std::vector<double>> nodes_temp, std::vector<std::vector<double>> weights_temp, int d, std::vector<int> Nl, std::vector<int> ids, double (*function_to_integrate)(std::vector<double> x, Args... rest), Args... rest);

void write_quadrature_points_to_file(std::ofstream& myfile, int iteration, std::vector<std::vector<double>> nodes_temp, int d, std::vector<int> Nl, std::vector<int> ids);
std::vector<double> van_der_corput_sequence(int p, int n, double epsilon);

template<typename... Args>
double integrate_with_sparse_grid(double (*function)(std::vector<double> x, Args... rest),int d, int l, std::vector<std::vector<double> > nodes,std::vector<std::vector<double> > weights,bool write_in_file, bool use_trap_rule,  Args... rest);

void monte_carlo_multivariate(std::vector<std::vector<double>>* nodes, std::vector<double>* weights, int l, int d, gsl_rng* r);

template<typename... Args>
double integrate_by_point_evaluation_multivariate(double (*function)(std::vector<double>* x, Args... rest),int n, int d, std::vector<std::vector<double>>* nodes, std::vector<double>* weights, Args... rest);



#endif /* multivariate_integration_hpp */
