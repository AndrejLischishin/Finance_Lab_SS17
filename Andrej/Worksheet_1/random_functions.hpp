//
//  reject_sample.hpp
//  Finanzpraktikum_W1
//
//  Created by WhoAmI on 22.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//

// enable doxygen processing for this header:
/** @file */
#ifndef random_functions_hpp
#define random_functions_hpp

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_randist.h>

double random_number01();
double random_number_01_GSL(gsl_rng* r);
void rejection_sampl_algo(int number_samples);

#endif /* random_functions_hpp */
