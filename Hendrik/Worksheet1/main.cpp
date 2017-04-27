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
	printf("Random Number: %lf\n", randomNumber01());
    //Print random number create by function randomNumber01GSL
	printf("Random Number GSL: %lf\n", randomNumber01GSL(r));

    //Open or create file "values.txt" to insert the values calculated by the rejection sampling algorithm    
    FILE *f = fopen("values.txt", "w");
    //Write 1000000 values into file
    for(int i=0; i<1000000; i++)
        fprintf(f, "%lf\n", rejectionSampling(r));
    //Close file "values.txt"
    fclose(f);

    //Task 4
    printf("Normal CDF at 0.5: %lf\n", normalCDF(0.5));
    //Open or create file "values2.txt" to insert the values calculated by the normalInverseCDF function    
    FILE *f2 = fopen("values2.txt", "w");
    for(int i=0; i<1000000; i++)
        fprintf(f2, "%lf\n", normalInverseCDF(randomNumber01GSL(r)));
    fclose(f2);

    //Task 6
    double *z = boxMuller(r);
    printf("Box Muller 1: %lf\n", z[0]);
    printf("Box Muller 2: %lf\n", z[1]);

    //Task 8
    int N = 1000000;
    double values[N];
    for(int i=0; i<N; i++)
        values[i] = rejectionSampling(r);
    
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

    printf("Sigma via algorithm: %lf\n", sigmaAlgorithm(values, N));
        
    printf("Mu: 5\n");

    //Task 9
	//Clean the files
	FILE *errorFile1 = fopen("error1.txt", "w");
    FILE *errorFile2 = fopen("error2.txt", "w");
    FILE *errorFile3 = fopen("error3.txt", "w");
    for(int k=1; k<=10000; k*=10)
    {
        double newValues[k];
        printf("Given Sigma: 0.1\n");
        for(int i=0; i<k; i++)
            newValues[i] = normalFromStandardNormalDistribution(r, 5., 0.1);
        printf("Sigma N=%i: %lf\n", k, sigmaAlgorithm(newValues, k));
        fprintf(errorFile1, "%i     %lf\n", k, fabs(sigmaAlgorithm(newValues, k)-0.1));

        printf("Given Sigma: 2\n");
        for(int i=0; i<k; i++)
            newValues[i] = normalFromStandardNormalDistribution(r, 5., 2);
        printf("Sigma N=%i: %lf\n", k, sigmaAlgorithm(newValues, k));
        fprintf(errorFile2, "%i     %lf\n", k, fabs(sigmaAlgorithm(newValues, k)-2.0));

        printf("Given Sigma: 10\n");
        for(int i=0; i<k; i++)
            newValues[i] = normalFromStandardNormalDistribution(r, 5., 10);
        printf("Sigma N=%i: %lf\n\n", k, sigmaAlgorithm(newValues, k));
        fprintf(errorFile3, "%i     %lf\n", k, fabs(sigmaAlgorithm(newValues, k)-10.0));
    }
    fclose(errorFile1);
    fclose(errorFile2);
    fclose(errorFile3);

    //Free allocated space for gsl
    gsl_rng_free(r);

	return 0;
}
