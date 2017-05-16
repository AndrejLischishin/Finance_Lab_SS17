#include <iostream>
#include "Sheet2.h"

int main(){

	// Task 9
	
	int N = 8;
	
	double* mc_val = new double[N];
	double* tr_val = new double[N];
	double* cc_val = new double[N];
	double* gl_val = new double[N];
	
	for(int i=1; i<N+2; i++){
		mc_val[i-1] = monte_carlo(i);
		tr_val[i-1] = trap_rule(i);
		cc_val[i-1] = clenshaw_curtis(i);
		gl_val[i-1] = gauss_legendre(i);
	}
	
	FILE *fp;
	fp = fopen("integration", "w");
	for(int i=1; i<N; i++){
	fprintf(fp,"%i  %f  %f  %f  %f\n", i, 
			mc_val[i-1], tr_val[i-1], cc_val[i-1], gl_val[i-1]);
		}
	fclose(fp);

	return 0;
}