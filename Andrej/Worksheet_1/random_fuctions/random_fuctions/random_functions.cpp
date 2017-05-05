//
//  random_functions.cpp
//  test
//
//  Created by WhoAmI on 28.04.17.
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
 * Rejection Sampling Algorithm
 *
 * The algorithm produces standard normal distributed value
 *
 * @param r a pointer to gsl_rng object
 * @return x returns double value, which is standard normal distributed
 *
 */

//rejection sampling algorithm
double rejection_sampl_algo(gsl_rng* r){
    
    /**
     * interval bounds \f$[a,b]\f$, s.t. \f$\int ^{b}_{a}p\left(x\right) dx = 1\f$,
     * p(x) density for a standard normal distribution
     */
    
    int a = -6, b = 6;
    //sigma = 1 beacuse standard normal distributed
    double sigma = 1;
    
    double x;
    double y ;
    
    //computing max_p(x),xє[a;b], max will be reached at x = 0, because normal distributed
    double max_px = gsl_ran_gaussian_pdf(0,sigma);
    double p_x;
    
    //until y>p(x), so until we are above the maximum of probability density function
    while (1) {
        //draws uniformly distributed x'є[a,b]
        x = gsl_ran_flat(r,a,b);
        //computes density value at point x
        p_x = gsl_ran_gaussian_pdf(x,sigma);
        //draws uniformly distributed yє[0,max_px]
        y = gsl_ran_flat(r,0,max_px);
        //check if the sampled point is under p(x), if so then return
        if(y<p_x||y==p_x){
            
            return x;
        }
    }
}

/**
 * Moro's algorithm is an approximation to the c.d.f. of the 
 * standard normal distribution with an accurancy of 8 digits
 *
 * @param x double value, point at which \f$p(x)\f$ will be calculated
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
 * Box Mueller Algorithm
 * @param mu It is mean for the normal distribution
 * @param sigma It is sigma for the normal distribution
 *
 * @return It returns normal distributed value
 *
 */

double mueller_box_algo(double mu, double sigma){
    
    const double epsilon = std::numeric_limits<double>::min();
    const double two_pi = 2.0 * 3.14159265358979323846;
    
    static double z0, z1;
    static bool generate;
    generate = !generate;
    
    
    if(!generate)
        return z1 * sigma + mu;
    
    double u1, u2;
    do
    {
        u1 = rand() * (1.0 / RAND_MAX);
        u2 = rand() * (1.0 / RAND_MAX);
    }
    while ( u1 <= epsilon );
    
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma + mu;
    
}

/**
 * naively computing mean mu and sigma of a collection
 * of normal distributed samples
 * @param N Number of given samples
 * @param sample Pointer to the vector of double valued
 * samples
 *
 * @return It returns sigma
 *
 */
double sigma_naive(int N, std::vector<double>* sample){
    
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




