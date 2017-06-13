#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include "integration_functions.hpp"
#include "simulation_functions.hpp"
#include "random_functions.hpp"


double testfunction(std::vector<double> x){
//	std::cout<<x[0]*x[0]+x[1]<<std::endl;
//	std::cout<<x[1]<<"x0"<<std::endl;
	return x[0]*x[0]+x[1];
}

double function_task13(std::vector<double>* x, double gamma, int d){
//	std::cout<<x[0]*x[0]+x[1]<<std::endl;
//	std::cout<<x[1]<<"x0"<<std::endl;
	double product = 1;
	for(int i=1; i<=d; i++){
		product = product*(1+gamma*exp((*x)[i-1]/2));
	}
	
	return product;
}

double function_task13_sparse(std::vector<double> x, double gamma, int d){
//	std::cout<<x[0]*x[0]+x[1]<<std::endl;
//	std::cout<<x[1]<<"x0"<<std::endl;
	double product = 1;
	for(int i=1; i<=d; i++){
		product = product*(1+gamma*exp(x[i-1]/2));
	}
	
	return product;
}

void monte_carlo(std::vector<std::vector<double>*>* nodes, std::vector<double>* weights, int l, int d, gsl_rng* r)
{
	int Nl = pow(2,l)-1;
	for(int i=1; i<=Nl; i++)
	{
		for(int j=1;j<=d;j++){
		(*nodes)[j-1]->push_back(random_number_01_GSL(r));
		}	
		(*weights)[i-1] = (1./(Nl));
	}
}

template<typename... Args>
double integrate_by_point_evaluation(double (*function_to_integrate)(std::vector<double>* x, Args... rest), std::vector<std::vector<double>*>* nodes, std::vector<double>* weights, int n, int d, Args... rest)
{
	double result = 0.0;
	for(int i=0; i<n; i++)
	{	
		result += (*weights)[i]*function_to_integrate((*nodes)[i], rest...);
	}

	return result;
}


void trap_rule_weights(std::vector<double>* weights, int l){
	int Nl = pow(2, l)-1;
	int Nk = pow(2, l-1)-1;
	for(int i=1; i<Nl-1; i++){
		 (*weights)[i] = (double) 1/(Nl+1);
		 if((i+1)%2==0){
		 	if(i==1 || i==Nl-2){
		 	(*weights)[i] = (*weights)[i]-(double) 3/(2*(Nk+1));
		 	}
		 	else (*weights)[i] = (*weights)[i]-(double) 1/(Nk+1);
		 }
	}
	(*weights)[0] = (double) 3/(2*(Nl+1));
	(*weights)[Nl-1] = (double) 3/(2*(Nl+1));
	
	if(l==1) (*weights)[0]=1;
	if(l==2) (*weights)[1]=-0.75;
}

void trap_rule_nodes(std::vector<double> *nodes, int l){
	int Nl = pow(2, l)-1;
	
	for(int i=1; i<=Nl; i++){
		 (*nodes)[i-1] = (double) i/(Nl+1);
	}
//	std::cout<<(*nodes)[0]<<"nodes "<<std::endl;
}

void clenshaw_curtis_nodes(std::vector<double>* nodes, int l)
{
	unsigned int Nl = pow(2,l)-1;

	for(unsigned int i=1; i<=Nl; i++)
	{
		(*nodes)[i-1]=(.5*(1.-cos((double)(M_PI*i/(Nl+1.)))));
		double sum = 0.;

	}
}


void clenshaw_curtis_weights(std::vector<double>* weights, int l){
	int Nl = pow(2, l)-1;
	int Nk = pow(2, l-1)-1;
	double sum = 0.;
	std::vector<double> weightstemp;

	for(unsigned int i=1; i<=Nl; i++)
	{
		sum = 0.;
		for(unsigned int j=1; j<=(Nl+1)/2; j++)
		{
			sum = sum + (double) (1./(2.*j-1.)) * sin( (double)(2.*j-1.) * M_PI * (i/(double)(Nl+1.)) );
		}
		(*weights)[i-1]=( (double)(2./(Nl+1.)) * sin( M_PI * ((double)i/(double)(Nl+1.)) ) * sum );
	}
	
	for(unsigned int i=1; i<=Nk; i++)
	{
		sum = 0.;

		for(unsigned int j=1; j<=(Nk+1)/2; j++)
		{
			sum = sum + (double) (1./(2.*j-1.)) * sin( (double)(2.*j-1.) * M_PI * (i/(double)(Nk+1.)) );
		}
		weightstemp.push_back( (double)(2./(Nk+1.)) * sin( M_PI * ((double)i/(double)(Nk+1.)) ) * sum );
	}
	
	for(int i=1; i<=Nl; i++){
		if((i)%2==0){
			(*weights)[i-1]=((*weights)[i-1]-weightstemp[(int)i/2-1]);
		//	(*weights)[i-1]=(*weights)[i-1]-;
		}
	//	std::cout<<(*weights)[i-1]<<" weightss "<<std::endl;
	} 
}


int enumeration(int *k, int d, int l, std::vector<int>* diag){
	int S_k=d;
	int II = 0;
	int counter=0;
	while(1<2){
	
/*	FILE *fp;
	fp = fopen("combinations", "a");
	
	for(int i = 0; i<pow(2,k[0])-1; i++){
		for(int j=0; j<pow(2,k[1])-1; j++){
			fprintf(fp,"%f  %f\n", nodes1[i], nodes2[j]);
		}
	}
	fclose(fp);
*/	
	for(int j=0; j<d; j++){
//		printf(" %i  enum ", k[j]);
//		diag[II+j] = k[j]-1;
		(*diag).push_back(k[j]-1);
	}
	II = II+d;
	counter++;

	// zur Kontrolle, damit man sehen kann, welche Kombinationen der Algorithmus erzeugt
	for(int i=0; i<d; i++){
//		printf("%i ", k[i]);
	}
//	printf("end \n"); 
	
	// Hier ist es wieder einfach der Algorithmus vom Arbeitsblatt
	
		for(int j=1; j<=d; j++){
			k[j-1]=k[j-1]+1;
			S_k = S_k +1;
			if(S_k > d+l-1){
				if(j==d){ return counter; }
				S_k = S_k-k[j-1]+1;
				k[j-1] = 1;
			}
			else break;
		}
		
	}
	
	return counter;
}


void loop(int* vec, int* klevel, int d, int* finalvec){
	int I;
	int count;
	int number=1;
	int II = 0;
	int last = d-1;
	
	while(last+1==d){
	
	for(int j=0; j<d; j++){
//		printf(" %i", vec[j]);
		finalvec[II+j] = vec[j]-1;
	}
	II = II+d;

//	printf(" end \n");
	number++;
		if(vec[last]<klevel[last]){
			vec[last]++;
			}
		else {
		for(int i=d-1; i>=0; i--){
			if(vec[i]<klevel[i]) {
			vec[i]++;
			break;
			}
			
			count = 0;
			for(I=0; I<d; I++){
			if(vec[I]==klevel[I]) { count++; 
				}
			}
			if(count==d) d=-2; 
			
			if(vec[i]==klevel[i]) vec[i]=1;
			}
		}
	}
}

void sparse_grid_nodes(int d, int product, int* allvec, std::vector<std::vector<double> > nodesv){
	FILE *fp;
	fp = fopen("stuetzstellen", "a");
		for(int i=0; i<product*d; i=i+d){
			for(int j=0; j<d; j++){
				fprintf(fp,"%f ", nodesv[j][allvec[i+j]]);
			}
		fprintf(fp,"\n");	
		}
	fclose(fp);	
}
	
void sparse_grid_weights(int d, int product, int* allvec, 
							double** weights, std::vector<std::vector<double> > weightsv){
	double product_w=1;
	FILE *fp;
	fp = fopen("gewichte", "a");
		for(int i=0; i<product*d; i=i+d){
			for(int j=0; j<d; j++){
				product_w=product_w*weights[j][allvec[i+j]];
				weightsv[j][allvec[i+j]]=weights[j][allvec[i+j]];
			}
			fprintf(fp,"%f \n ", product_w);
		fprintf(fp,"\n");	
		}
	fclose(fp);	
}
	
	
template<typename... Args>
double integrate_with_sparse_grid(double (*function_to_integrate)(std::vector<double> x, Args... rest), 
			int d, int l, std::vector<std::vector<double> > nodes, 
			std::vector<std::vector<double> > weights, bool write_in_file, bool use_trap_rule, Args... rest){
	if(write_in_file==true) {
		FILE *fp;
	fp = fopen("stuetzstellen", "w");
		fprintf(fp,"");	
		
	fclose(fp);
	
	}
	int maxlevel = (int)pow(2,l)-1;
	int* klevel = new int[d];
	int* k = new int[d];
	int* k1 = new int[d];
	int* vec = new int[d];
	int sum = 1;
	

//	int sz = (l*(l+1)*0.5)*d;
//	sz = 6*d;
//	int* diag = new int[sz];
	std::vector<int> diag;
	int K;
	double final_sum;
	int J;
	double product_w=1;
	int product = 1;
	double final_value=0;
	
	for(int i=0; i<d; i++){
		vec[i]=1;
		k1[i]=1;
		sum = sum * maxlevel;
	}
	
	printf("%i maxlevel \n", maxlevel);
	int* allvec = new int[sum*d];
	
	int sz = enumeration(k1,d,l,&diag);
//	std::cout<<sz<<std::endl;
	for(K=0; K<sz; K++)
//	while(diag[K+d+1]!=-1)
	{
	// von diag auf k uebertragen
		for(int i=0; i<d; i++){
			k[i]=diag[i+K*d]+1;
	//		printf("%i ki \n", k[i]);
		}
	
	
	// klevel
		for(int i=0; i<d; i++){
			klevel[i]=pow(2,k[i])-1;
		//	printf("%i klevel \n", klevel[i]);
		}
	
	
	
	// vec auf 1 setzen
		for(int i=0; i<d; i++){
			vec[i]=1;
		}

		// Tensorprodukt 
		loop(vec, klevel, d, allvec);
	
		// produkt berechnen
		product = 1;
		for(int i=0; i<d; i++){
			product = product * klevel[i];
		}
		
		if(use_trap_rule== true){
		for(int i=0; i<d; i++){
		trap_rule_weights(&weights[i], k[i]);
		trap_rule_nodes(&nodes[i], k[i]);
		}
		}
		else{
		for(int i=0; i<d; i++){
		clenshaw_curtis_weights(&weights[i], k[i]);
		clenshaw_curtis_nodes(&nodes[i], k[i]);
		}
		}
		
//		std::cout<<nodes[1][2]<<" nodes "<<std::endl;
		for(int i=0; i<product*d; i=i+d){
			for(int j=0; j<d; j++){
				product_w = product_w*weights[j][allvec[i+j]];
				J =j;
			//	std::cout<<weights[j][allvec[i+j]]<<std::endl;
			}
		//nodes[J][allvec[i]]*nodes[J][allvec[i]]+nodes[J][allvec[i+1]]
			std::vector<double> point;
			point.clear();
			for(int h=0; h<d; h++){
				point.push_back(nodes[h][allvec[i+h]]);
		//		std::cout<<point[h]<<" point "<<std::endl;
			}
			final_value = final_value + product_w*(function_to_integrate(point, rest...));
			product_w = 1;
	//		std::cout<<product<<std::endl;
	//		std::cout<<final_value<<std::endl;
	
	if(write_in_file==true) sparse_grid_nodes(d,product,allvec,nodes);
			
		}
			
	}
	free(allvec);
	free(k1);
	free(k);
	free(vec);	
	
	return final_value;		
	
}



	
int main(){


	int d=7;
	int l=4;
	int maxlevel = pow(2,l)-1;
	
	
/*	
	double** nodes = new double*[d];
	double** weights = new double*[d];
	for(int i=0; i<d; i++) {
		weights[i] =  new double[maxlevel];
		nodes[i] = new double[maxlevel];
	}	
*/	
	std::vector<std::vector<double> > nodesv(d);
	for ( int i = 0 ; i < d ; i++ ){
   		nodesv[i].resize(maxlevel);
	}
	
	std::vector<std::vector<double> > weightsv(d);
	for ( int i = 0 ; i < d ; i++ ){
   		weightsv[i].resize(maxlevel);
	}


//	sparse_grid_nodes(d, product, allvec, nodes, nodesv);
//	sparse_grid_weights(d, product, allvec, weights, weightsv);
	
	double final_value = integrate_with_sparse_grid(function_task13_sparse, 
			d, l, nodesv, weightsv, false, false, 0.1, d);
	
	
	
	printf("%f final sum\n", final_value);
	

}







