#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "Sheet1.h"

// Calculation of weights and abscissae

void trap_rule_nw(double* nodes, double* weights, int l){
	int Nl = pow(2, l)-1;
	for(int i=1; i<=Nl; i++){
		 nodes[i-1] = (double) i/(Nl+1);
		 weights[i-1] = (double) 1/(Nl+1);
	}
	weights[0] = (double) 3/(2*(Nl+1));
	weights[Nl] = (double) 3/(2*(Nl+1));
}


void clenshaw_curtis_nw(double* nodes, double* weights, int l){
	int Nl = pow(2, l)-1;
	for(int i=1; i<Nl; i++){
	
		 nodes[i-1] = .5*(1-cos((double)(M_PI*i/(Nl+1))*M_PI/180));
		 
		 	double sum = 0;
			for(int j=1; j<=(Nl+1)/2; j++){
			sum = sum + (double) 1/(2*j-1)*sin((double)(2*j-1)*M_PI*i/(Nl+1)*M_PI/180);
		}
		 weights[i-1] = (double)2/(Nl+1)*sin((double)M_PI*i/(Nl+1)*M_PI/180)*sum;
	}
	
}


void gauss_legendre_nw(double *nodes, double* weights, size_t l){
	
	double* wi = new double[1];
	double* xi = new double[1];
	
	double a=0; 
	double b=1;
	
	gsl_integration_glfixed_table* table;	
	table =	gsl_integration_glfixed_table_alloc(l);

	for(int j=0; j<l; j++){
	gsl_integration_glfixed_point (a, b, j, xi, wi,table);
	nodes[j]=xi[0];
	weights[j]=wi[0];
	}
}

// calculation of the integral

double monte_carlo(int N){
	
	int Nl = pow(2, N)-1;
	double sum = 0;
	
	gsl_rng* r;
	r = gsl_rng_alloc(gsl_rng_mt19937);
	
	for(int i=1; i<=Nl; i++){
		sum = sum + 1 + exp(.5*randomNumber01GSL(r));
	}
	return sum/(Nl);
}


double trap_rule(int N){
	int Nl = pow(2, N)-1;
	double* nodes = new double[Nl];
	double* weights = new double[Nl];
	
	trap_rule_nw(nodes, weights, N);
	
	double sum = 0;
	for(int i=0; i<Nl; i++){
		sum = sum + weights[i]*(1+exp(0.5*nodes[i]));
	}
	return sum;
}


double clenshaw_curtis(int N){
	int Nl = pow(2, N)-1;
	double* nodes = new double[Nl];
	double* weights = new double[Nl];
	
	clenshaw_curtis_nw(nodes, weights, N);
	
	double sum = 0;
	for(int i=0; i<Nl; i++){
		sum = sum + weights[i]*(1+exp(0.5*nodes[i]));
	}
	return sum;
}

	
double gauss_legendre(int N){
	size_t Nl = pow(2, N)-1;
	
	double* nodes = new double[Nl];
	double* weights = new double[Nl];
	
	gauss_legendre_nw(nodes, weights, Nl);
	
	double sum = 0;
	for(int i=0; i<Nl; i++){
		sum = sum + weights[i]*(1+exp(0.5*nodes[i]));
	}
	return sum;
}




















