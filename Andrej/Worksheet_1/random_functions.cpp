//
//  reject_sample.cpp
//  Finanzpraktikum_W1
//
//  Created by WhoAmI on 22.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//

/** @file */
#include "random_functions.hpp"

/**
 * Draws a random number bewtween \f$[0,1]\f$ via rand.
 *
 * @return The drawn random number.
 */
double random_number01()
{
	return (double) rand() / RAND_MAX; // *
}

/**
 * Draws a random number in \f$[0,1]\f$ via the gsl.
 *
 * @param r A pointer to the random number generator which is used.
 *
 * @return The drawn random number.
 */
double random_number_01_GSL(gsl_rng* r)
{
	return gsl_rng_uniform(r);
}

// enable doxygen processing for this header:
/** @file */

/**
 * Main function for random number evaluation.
 * A first random number \f$x_1\f$ is drawn via rand.
 * A second random number \f$x_2\f$ is drawn via gsl_rng_uniform.
 *
 * Furthermore, an array of 10 doubles is allocated. However, someone
 * seems to have forgotten to free the allocated space again in the end...
 *
 * @param argc The number of arguments provided.
 * @param argv An array of arguments (argv[0] is the name of the
 * executable).
 *
 * @return If everything worked fine, 0 is returned.
 */


/**
 * rejection sampling algorithm
 * this function produces specified by the input number of standard normal
 * distributed values, which will be automatically written to the
 * "rejection_sampl.txt" file
 * @param number_samples an integer argument, specifies number of velues, default = 1
 * 
 *
 */


void rejection_sampl_algo(int number_samples = 1){

  /**
   * interval bounds \f$[a,b]\f$, s.t. \f$\int ^{b}_{a}p\left(x\right) dx = 1\f$,
   * p(x) density for a standard normal distribution
   */
    int a = -6, b = 6;
    double sigma = 1;
    double x;
    double y;
    gsl_rng* r;

    //computing max_p(x),xє[a;b], max will be reached at x = 0, because normal distributed
    double p_x;
    //memory allocation
    r = gsl_rng_alloc(gsl_rng_mt19937);

    //integer for iterating
    int ittr = 1;
    std::ofstream myfile;
    //opens a file and checks if it was successfully
    //when file will be opened previous content will be deleted
    myfile.open("rejection_sampl.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    while (ittr<=number_samples){
      //until y>p(x), so until we are above the maximum of probability density function
      while (1) {
          //generates uniformly(flat) distributed value
          x = gsl_ran_flat(r,a,b);
          //computes density value at point x, with sigma = 1
          p_x = gsl_ran_gaussian_pdf(x,sigma);
          //check if the sampled point is under p(x), if so then break
          if(y<p_x||y==p_x){
              //writes values into a file
              myfile<<x<<std::endl;
              break;
            }
          ittr++;
      }
    }

    //closes file
    myfile.close();
    //frees memory
    gsl_rng_free(r);
  }
