//
//  random_functions.hpp
//  test
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//

#ifndef random_functions_hpp
#define random_functions_hpp

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#define A0 0.398942270991
#define A1 0.020133760596
#define A2 0.002946756074
#define B1 0.217134277847
#define B2 0.018576112465
#define B3 0.000643163695
#define C0 1.398247031184
#define C1 -0.360040248231
#define C2 0.022719786588
#define D0 1.460954518699
#define D1 -0.305459640162
#define D2 0.038611796258
#define D3 -0.003787400686
#define E0 2.50662823884
#define E1 -18.61500062529
#define E2 41.39119773534
#define E3 -25.44106049637
#define F0 -8.47351093090
#define F1 23.08336743743
#define F2 -21.06224101826
#define F3 3.13082909833
#define G0 0.3374754822726147
#define G1 0.9761690190917186
#define G2 0.1607979714918209
#define G3 0.0276438810333863
#define G4 0.0038405729373609
#define G5 0.0003951896511919
#define G6 0.0000321767881768
#define G7 0.0000002888167364
#define G8 0.0000003960315187

double random_number01();
double random_number_01_GSL(gsl_rng* r);
double rejection_sampl_algo(gsl_rng* r);
double normal_cdf(double x);
double normal_inverse_cdf(double x);
std::vector<double>* box_muller_algo(gsl_rng* r);
double sigma_naive(std::vector<double>* sample, int N);
double sigma_algorithm(std::vector<double>* sample, int N);
std::vector<double>* brownian_motion(gsl_rng* r, double T, double delta_t, std::vector<double>* w,double s0, double mu, double sigma);
std::vector<double>* wiener_process(gsl_rng* r, double T, double delta_t);


#endif /* random_functions_hpp */
