#include<iostream>
#include"Sheet1.h"

int main() {

	gsl_rng* r;
	r = gsl_rng_alloc(gsl_rng_mt19937);
	
	// Task 2: Rejection Sampling Algorithm 
	
	double a = -2;
	double b = 2;
	double max = 1./sqrt(2*M_PI);
	
	FILE *fp1;
	fp1 = fopen("values.txt", "w");
	for(int i=1; i<10000; i++){
	fprintf(fp1,"%f \n", rej_sam_alg(a, b, max, r));
	}
	fclose(fp1);
	
	// Task 6 Box Muller Method

	int n=1000;
	double* z_1 = (double*)malloc(n*sizeof(double));
	double* z_2 = (double*)malloc(n*sizeof(double));
	
	
	box_muller(z_1, z_2, n, r);
	
	FILE *fp2;
	fp2 = fopen("box_muller", "w");
	for(int i=1; i<n; i++){
	fprintf(fp2,"%f %f\n", z_1[i], z_2[i]);
	}
	fclose(fp2);

	// Task 9 Parameter Estimation
	
	int N=10000;
	
	double* z_1N = (double*)malloc(N*sizeof(double));
	double* z_2N = (double*)malloc(N*sizeof(double));

	double ex = 2.0;
	
	double sigma = 0.1;

	FILE *fp3;
	fp3 = fopen("variance_0.1", "w");
	for(int i=10; i<N; i++){
	box_muller_ex_var(z_1N, z_2N, i, r, ex, sigma);
	fprintf(fp3,"%f \n", fabs(var_est(z_1, i)-sigma));
	}
	
	fclose(fp3);
	
	sigma = 1.0;
	
	FILE *fp4;
	fp4 = fopen("variance_1.0", "w");
	for(int i=10; i<N; i++){
	box_muller_ex_var(z_1N, z_2N, i, r, ex, sigma);
	fprintf(fp4,"%f \n", fabs(var_est(z_1, i)-sigma));
	}
	
	fclose(fp4);
	
	sigma = 10.0;
	
	FILE *fp5;
	fp5 = fopen("variance_10.0", "w");
	for(int i=10; i<N; i++){
	box_muller_ex_var(z_1N, z_2N, i, r, ex, sigma);
	fprintf(fp5,"%f \n", fabs(var_est(z_1, i)-sigma));
	}
	
	fclose(fp5);

	// Task 10: Wiener Process
	
	double T = 2;
	double delta_t = 0.01;
	
	double M = (double) T/delta_t;
	double* w = (double*)malloc(M*sizeof(double));
	
	// Simulation fuer obenstehende Werte
	
	wiener_process(w, 0.1, 0.2, delta_t, T, r);

	FILE *fp6;
	fp6 = fopen("wiener_process1", "w");
	for(int i=0; i<M; i++){	
	fprintf(fp6,"%f \n", w[i]);
	}
	fclose(fp6);
	
	wiener_process(w, 0.1, 0.2, delta_t, T, r);
	
	FILE *fp7;
	fp7 = fopen("wiener_process2", "w");
	for(int i=0; i<M; i++){	
	fprintf(fp7,"%f \n", w[i]);
	}
	fclose(fp7);
	
	wiener_process(w, 0.1, 0.2, delta_t, T, r);
	
	FILE *fp8;
	fp8 = fopen("wiener_process3", "w");
	for(int i=0; i<M; i++){	
	fprintf(fp8,"%f \n", w[i]);
	}
	fclose(fp8);
	

}