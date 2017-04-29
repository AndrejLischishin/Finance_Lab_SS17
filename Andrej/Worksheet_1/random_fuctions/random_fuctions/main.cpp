//
//  main.cpp
//  test
//
//  Created by WhoAmI on 28.04.17.
//  Copyright © 2017 Andrei. All rights reserved.
//

#include "random_functions.hpp"

int main(int argc, char* argv[]){
    
    gsl_rng* r;
    
    // seeding
    unsigned long seed = std::time(NULL);
    //memory allocation
    r = gsl_rng_alloc(gsl_rng_mt19937);

    ///////////////////////////////////////////
    ///////////////////Task_1//////////////////
    ///////////////////////////////////////////
    
    gsl_rng_set(r, seed);
    
    printf( "%lf\n", random_number01());
    
    printf( "%lf\n", random_number_01_GSL(r));
    
    ///////////////////////////////////////////
    ///////////////////Task_2//////////////////
    ///////////////////////////////////////////
    
    int number_samples = 1000000;
    
    std::ofstream myfile;
    
    //opens a file and checks if it was successfully
    //when file will be opened previous content will be deleted
    
    myfile.open("rejection_sampl.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    
    //for iteration
    int ittr = 1;
    
    while (ittr<=number_samples) {
        myfile<<rejection_sampl_algo(r)<<std::endl;
        ittr++;
    }
    //closes file
    myfile.close();
    
    
    
    ///////////////////////////////////////////
    ///////////////////Task_3_matlab///////////
    ///////////////////////////////////////////
    
    ///////////////////////////////////////////
    ///////////////////Task_4_5////////////////
    ///////////////////////////////////////////
    
    /**
     * Beasley_Springer_Moro Algorithm
     * it computes inverse c.d.f. for standard normal distribution 
     * at x in order to become a standard normal distributed value
     *
     */
    normal_inverse_cdf(gsl_ran_flat(r, 0, 1));
    
    //frees memory
    gsl_rng_free(r);
    
    ///////////////////////////////////////////
    ///////////////////Task_6_7////////////////
    ///////////////////////////////////////////
    
    double mu = 0;
    double sigma = 1;
    
    //write 1000 samples for 2D plot in the file
    int num_sampl = 1000;
    srand((unsigned)std::time(NULL));
    
    //opens a file and checks if it was successfully
    //when file will be opened previous content will be deleted
    myfile.open("mueller_box.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    
    for (int i = 0; i<num_sampl; i++) {
        myfile<<mueller_box_algo(mu, sigma)<<std::endl;
    }
        
    //closes file
    myfile.close();

    ///////////////////////////////////////////
    ///////////////////Task_8_9////////////////
    ///////////////////////////////////////////
    
    
    int N_max = 10000000;
    
    //Three different values for sigma a
    int num_of_sigmas = 3;
    double sigma_s[] = {0.1,1,10};
    double mean = 0;
    
    std::string files[] = {std::to_string((int)sigma_s[0]),std::to_string((int)sigma_s[1]),std::to_string((int)sigma_s[2])};
    std::vector<double> samples;
    
    double sigma_approx = 0.0;
    std::vector<std::vector<double>> sigma_err(num_of_sigmas);
    
    srand((unsigned)std::time(NULL));
    
    
    for (int j = 0; j<num_of_sigmas; j++) {
        
        //opens a file and checks if it was successfully
        //when file will be opened previous content will be deleted
        myfile.open("sigma_err_"+files[j]+".txt",std::ios::trunc);
        //if (!myfile.is_open()) {
          //  std::cout<<"Error opening the file"<<std::endl;
        //}
        
        for (int N = 10,k = 0; N<=N_max; N = 10*N, k++) {
            
            for (int i = 0; i<N; i++) {
                
                samples.push_back(mueller_box_algo(mean, sigma_s[j]));
                
            }
            
            sigma_approx = sigma_algorithm(&samples, N);
            sigma_err[j].push_back(fabs(sigma_approx - sigma_s[j]));
            myfile<<sigma_err[j][k]<<" "<<N<<std::endl;
            samples.clear();
        }
        
        //closes file
        myfile.close();
        
    }
    

    return 0;
    
}










