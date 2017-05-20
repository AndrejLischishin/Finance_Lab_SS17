//
//  simulation_functions.hpp
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//

#ifndef simulation_functions_hpp
#define simulation_functions_hpp

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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

std::vector<double>* brownian_motion(gsl_rng* r, double T, double delta_t, std::vector<double>* w,double s0, double mu, double sigma);
std::vector<double>* wiener_process(gsl_rng* r, double T, double delta_t);

#endif /* simultaion_functions_hpp */
