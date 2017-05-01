#include "random_functions.hpp"

// enable doxygen processing for this header:
/** @file */

/** 
 * Executes all the Tasks from the first worksheet.
 *
 * @param argc The number of arguments provided.
 * @param argv An array of arguments (argv[0] is the name of the
 * executable).
 *
 * @return If everything worked fine, 0 is returned.
 */
int main(int argc, char* argv[]) 
{
	gsl_rng* r;
	r = gsl_rng_alloc(gsl_rng_mt19937);
    struct timeval tv;
    gettimeofday(&tv,0);
    unsigned long seed = tv.tv_sec;
    gsl_rng_set(r, seed);
    //Print random number created by function randomNumber01
	printf("Random Number: %lf\n", random_number_01());
    //Print random number create by function randomNumber01GSL
	printf("Random Number GSL: %lf\n", random_number_01_GSL(r));

    //Open or create file "values.txt" to insert the values calculated by the rejection sampling algorithm    
    FILE *f = fopen("values.txt", "w");
    //Write 1000000 values into file
    for(int i=0; i<1000000; i++)
        fprintf(f, "%lf\n", rejection_sampling(r));
    //Close file "values.txt"
    fclose(f);

    //Task 4
    printf("Normal CDF at 0.5: %lf\n", normal_CDF(0.5));
    //Open or create file "values2.txt" to insert the values calculated by the normalInverseCDF function    
    FILE *f2 = fopen("values2.txt", "w");
    for(int i=0; i<1000000; i++)
        fprintf(f2, "%lf\n", normal_inverse_CDF(random_number_01_GSL(r)));
    fclose(f2);

    //Task 6
    double *z = box_muller(r);
    printf("Box Muller 1: %lf\n", z[0]);
    printf("Box Muller 2: %lf\n", z[1]);

    //Task 8
    int N = 1000000;
    double values[N];
    for(int i=0; i<N; i++)
        values[i] = rejection_sampling(r);
    
    double mu = 0.;
    for(int i=0; i<N; i++)
        mu += values[i];
    mu = mu/N;
    printf("Mu: %lf\n", mu);

    double sigma = 0.;
    for(int i=0; i<N; i++)
        sigma += pow((values[i]-mu), 2.);
    sigma = sigma/(N-1);
    printf("Sigma naively: %lf\n", sigma);

    printf("Sigma via algorithm: %lf\n", sigma_algorithm(values, N));
        
    printf("Mu: 5\n");

    //Task 9
	//Clean the files
	FILE *error_file_1 = fopen("error1.txt", "w");
    FILE *error_file_2 = fopen("error2.txt", "w");
    FILE *error_file_3 = fopen("error3.txt", "w");
    for(int k=1; k<=10000; k*=10)
    {
        double new_values[k];
        printf("Given Sigma: 0.1\n");
        for(int i=0; i<k; i++)
            new_values[i] = normal_from_standard_normal_distribution(r, 5., 0.1);
        printf("Sigma N=%i: %lf\n", k, sigma_algorithm(new_values, k));
        fprintf(error_file_1, "%i     %lf\n", k, fabs(sigma_algorithm(new_values, k)-0.1));

        printf("Given Sigma: 2\n");
        for(int i=0; i<k; i++)
            new_values[i] = normal_from_standard_normal_distribution(r, 5., 2);
        printf("Sigma N=%i: %lf\n", k, sigma_algorithm(new_values, k));
        fprintf(error_file_2, "%i     %lf\n", k, fabs(sigma_algorithm(new_values, k)-2.0));

        printf("Given Sigma: 10\n");
        for(int i=0; i<k; i++)
            new_values[i] = normal_from_standard_normal_distribution(r, 5., 10);
        printf("Sigma N=%i: %lf\n\n", k, sigma_algorithm(new_values, k));
        fprintf(error_file_3, "%i     %lf\n", k, fabs(sigma_algorithm(new_values, k)-10.0));
    }
    fclose(error_file_1);
    fclose(error_file_2);
    fclose(error_file_3);

	//Task 10
	double s0 = 10.0;
	mu = 0.1;
	sigma = 0.2;
	double T = 2;
	FILE *wiener_process_1 = fopen("wiener_process_1.txt", "w");
	FILE *wiener_process_2 = fopen("wiener_process_2.txt", "w");
	double delta_t = 0.5;
	int M = (int)(T/delta_t);
	double *w_1 = wiener_process(r, T, delta_t);
	double *w_2 = wiener_process(r, T, delta_t);
	double *w_3 = wiener_process(r, T, delta_t);
	for(int i=0; i<=M; i++)
		fprintf(wiener_process_1, "%lf     %lf		%lf		%lf\n", i*delta_t, w_1[i], w_2[i], w_3[i]);
	delta_t = 0.01;
	M = (int)(T/delta_t);
	double *w_4 = wiener_process(r, T, delta_t);	
	double *w_5 = wiener_process(r, T, delta_t);
	double *w_6 = wiener_process(r, T, delta_t);
	for(int i=0; i<=M; i++)
		fprintf(wiener_process_2, "%lf     %lf		%lf		%lf\n", i*delta_t, w_4[i], w_5[i], w_6[i]);
	fclose(wiener_process_1);	
	fclose(wiener_process_2);
	
	FILE *brownian_motion_1 = fopen("brownian_motion_1.txt", "w");
	FILE *brownian_motion_2 = fopen("brownian_motion_2.txt", "w");
	delta_t = 0.5;
	M = (int)(T/delta_t);
	double *s_1 = brownian_motion(r, T, delta_t, w_1, s0, mu, sigma);
	double *s_2 = brownian_motion(r, T, delta_t, w_2, s0, mu, sigma);
	double *s_3 = brownian_motion(r, T, delta_t, w_3, s0, mu, sigma);
	for(int i=0; i<=M; i++)
		fprintf(brownian_motion_1, "%lf     %lf		%lf		%lf\n", i*delta_t, s_1[i], s_2[i], s_3[i]);
	delta_t = 0.01;
	M = (int)(T/delta_t);
	double *s_4 = brownian_motion(r, T, delta_t, w_4, s0, mu, sigma);	
	double *s_5 = brownian_motion(r, T, delta_t, w_5, s0, mu, sigma);
	double *s_6 = brownian_motion(r, T, delta_t, w_6, s0, mu, sigma);
	for(int i=0; i<=M; i++)
		fprintf(brownian_motion_2, "%lf     %lf		%lf		%lf\n", i*delta_t, s_4[i], s_5[i], s_6[i]);
	fclose(brownian_motion_1);	
	fclose(brownian_motion_2);


    //Free allocated space for gsl
    gsl_rng_free(r);

	return 0;
}
