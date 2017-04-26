#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_randist.h>
#include<iostream>

// enable doxygen processing for this header:
///**
 //@file
 //*/

int main(int argc, char* argv[]){


  std::cout<<"rand x = "<<(double)rand() / RAND_MAX<<std::endl;

/**
 * (double) cust an integer into a double, so that we have a
 * double division, if we dont use it the output will be
 * always 0, because rand() / RAND_MAX is between 0 and 1
 */

  gsl_rng* r;
  r = gsl_rng_alloc(gsl_rng_mt19937);
  std::cout<<"gsl_uniform x = "<<gsl_rng_uniform(r)<<std::endl;
  gsl_rng_free(r);
 /**
  * Yes GSL has function which simulates normally distributed random variables
  * You need to include <gsl/gsl_randist.h> and <gsl/gsl_cdf.h> in order to use different functions
  * for Gaussian Distribution
  * mean 0, standard diviation 2
  */

  std::cout<<"gsl_gaussian x = "<<gsl_ran_gaussian(r,2)<<std::endl;

  return 0;

}
