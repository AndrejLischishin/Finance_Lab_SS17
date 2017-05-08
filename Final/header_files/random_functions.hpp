//
//  random_functions.hpp
//  test
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//

#ifndef random_functions_hpp
#define random_functions_hpp

/** @file */
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


//
//  random_functions.cpp
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//


/**
 * Draws a random number bewtween \f$[0,1]\f$ via rand.
 *
 * @return The drawn random number.
 */
double random_number_01()
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

/**
 * Rejection Sampling Algorithm
 *
 * The algorithm produces a standard normal distributed value.
 *
 * @param r A pointer to gsl_rng object
 * @return x Returns double value, which is standard normal distributed
 *
 */

//rejection sampling algorithm
double rejection_sampl_algo(gsl_rng* r)
{    
    /**
     * interval bounds \f$[a,b]\f$, s.t. \f$\int ^{b}_{a}p\left(x\right) dx = 1\f$,
     * p(x) density for a standard normal distribution
     */
    
    int a = -6, b = 6;
    //sigma = 1 beacuse standard normal distributed
    double sigma = 1;
    
    double x;
    double y;
    
    //computing max_p(x),xє[a;b], max will be reached at x = 0, because normal distributed
    double max_px = gsl_ran_gaussian_pdf(0,sigma);
    double p_x;
    
    //until y>p(x), so until we are above the maximum of probability density function
    do{
        //draws uniformly distributed x'є[a,b]
        x = gsl_ran_flat(r,a,b);
        //computes density value at point x
        p_x = gsl_ran_gaussian_pdf(x,sigma);
        //draws uniformly distributed yє[0,max_px]
        y = gsl_ran_flat(r,0,max_px);
        
        //check if the sampled point is under p(x), if so then return
    }while(y<=p_x);
    
    return x;
    
}

/**
 * Moro's algorithm is an approximation to the c.d.f. of the 
 * standard normal distribution with an accurancy of 8 digits.
 *
 * @param x Double value, point at which \f$p(x)\f$ will be calculated
 *
 * @return If everything worked fine returns \f$p(x)\f$
 */
double normal_cdf(double x){
    
    double x2;
    if(x<0.0)
        return 1.0-normal_cdf(-x);
    if(x<=1.87)
    {
        x2 = x*x;
        return 0.5+x*(A0+(A1+A2*x2)*x2)/(1.0+(B1+(B2+B3*x2)*x2)*x2);
    }
    else if(x<6.0)
    {
        return 1.0-pow((C0+(C1+C2*x)*x)/(D0+(D1+(D2+D3*x)*x)*x), 16.0);
    }
    return 1.0;
}


/**
 * Calculates the inverse CDF of the standard normal distribution for a parameter x.
 *
 * @param x The parameter for the inverse CDF of the standard normal distribution.
 *
 * @return The value of the inverse CDF at x.
 */

double normal_inverse_cdf(double x){
    
    double p = x-0.5;
    double r;
    if(fabs(p)<0.42)
    {
        r = pow(p, 2.);
        return p*(((E3*r+E2)*r+E1)*r+E0)/((((F3*r+F2)*r+F1)*r+F0)*r+1.0);
    }
    else
    {
        if(p<0)
            r = x;
        else
            r = 1-x;
        
        r = log(-log(r));
        r = G0+r*(G1+r*(G2+r*(G3+r*(G4+r*(G5+r*(G6+r*(G7+r*G8)))))));
        
        if(p<0)
            return -r;
        else
            return r;
    }
}

/**
 * Box Muller Algorithm.
 *
 * @param mu It is mean for the normal distribution
 * @param sigma It is sigma for the normal distribution
 *
 * @return It returns pointer to the vector with 2 normal distributed values
 *
 */

std::vector<double>* box_muller_algo(gsl_rng *r)
{        
    std::vector<double>* z;
    z = new std::vector<double>;
    
    
    double z0, z1;
    double u1 = random_number_01_GSL(r);
    double u2 = random_number_01_GSL(r);
    
    z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(2 * M_PI * u2);
    
    z->push_back(z0);
    z->push_back(z1);
    
    return z;
    
}

/**
 * naively computing mean mu and sigma of a collection
 * of normal distributed samples
 * @param N Number of given samples
 * @param sample Pointer to the vector of double valued
 * normal distributed samples
 *
 * @return It returns sigma
 *
 */
double sigma_naive(std::vector<double>* sample, int N)
{    
    double sigma = 0.0;
    double mu = 0.0;   
    
    for (int i = 0; i<N; i++) {
        mu += (*sample)[i];
    }
    
    mu = mu/(double)N;
    
    for (int i = 0; i<N; i++) {
        sigma += ((*sample)[i]-mu)*((*sample)[i]-mu);
    }
    
    sigma = sqrt((1./(double)(N-1))*sigma);
    
    return sigma;
}


/**
 * Calculates the variance for N given values.
 *
 * @param sample Samples to calculate the variance of.
 * @param N Number of samples.
 *
 * @return The calculated variance of the samples.
 */
double sigma_algorithm(std::vector<double>* sample, int N)
{
    double alpha = (*sample)[0];
    double beta = 0.;
    int i;
    double sigma = 1.;
    for(i=1; i<N; i++)
    {
        double gamma = (*sample)[i]-alpha;
        alpha = alpha+gamma/(i+1);
        beta = beta+gamma*gamma*i/(i+1);
    }
    
    sigma = sqrt(beta/i);
    
    return sigma;
}


/**
 * Simulates a wiener process.
 *
 * @param r Pointer to the gsl_rng object for generating standard normal distributed numbers
 * @param T Time period of simulated process
 * @param delta_t Step of discretisation
 *
 * @return Pointer to the vector of values at discretisation points
 */

std::vector<double>* wiener_process(gsl_rng* r, double T, double delta_t)
{
    int M = (int)(T/delta_t);
    std::vector<double> *w = new std::vector<double>(M+1);
    (*w)[0] = 0.0;

    for(int i=0; i<M; i++)
    {
        (*w)[i+1] = (*w)[i]+sqrt(delta_t)*gsl_ran_ugaussian(r);
    }
    
    return w;
}

/**
 * Simulates brownian motion path for the given values of wiener process.
 *
 * @param r Pointer to the gsl_rng object for generating standard normal distributed numbers
 * @param T Time period of simulated process
 * @param delta_t Step of discretisation
 * @param w Pointer to the vector with values of wiener process ar discretisation points
 * @param s0 Value of brownian_motion at time = 0
 * @param mu Drift
 * @param sigma Volatility
 *
 * @return Pointer to the vector of values at discretisation points
 *
 */
    
std::vector<double>* brownian_motion(gsl_rng* r, double T, double delta_t, std::vector<double>* w, double s0, double mu, double sigma)
{
    int M = (int)(T/delta_t);
    std::vector<double> *s = new std::vector<double>(M+1);
    (*s)[0] = s0;
    for(int i=0; i<M; i++)
    {
        (*s)[i+1] = s0*exp((mu-0.5*pow(sigma, 2.0))*(i+1)*delta_t+sigma*(*w)[i+1]);
    }
    
    return s;
}

#endif /* random_functions_hpp */
