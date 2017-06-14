
#include "multivariate_integration.hpp"

double testfunction(std::vector<double> x){
//	std::cout<<x[0]*x[0]+x[1]<<std::endl;
//	std::cout<<x[1]<<"x0"<<std::endl;
	return x[0]*x[0]+x[1];
}


void trap_rule_s(std::vector<double>* weights, int l){
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

void trap_rulen(std::vector<double> *nodes, int l){
	int Nl = pow(2, l)-1;

	for(int i=1; i<=Nl; i++){
		 (*nodes)[i-1] = (double) i/(Nl+1);
	}
//	std::cout<<(*nodes)[0]<<"nodes "<<std::endl;
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

void sparse_grid_nodes(int d, int product, int* allvec,
							double** nodes, std::vector<std::vector<double> > nodesv){
	FILE *fp;
	fp = fopen("stuetzstellen", "a");
		for(int i=0; i<product*d; i=i+d){
			for(int j=0; j<d; j++){
				fprintf(fp,"%f ", nodes[j][allvec[i+j]]);
				nodesv[j][allvec[i+j]]=nodes[j][allvec[i+j]];
			}
		fprintf(fp,"\n");
		}
	printf("\n");
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
double integrate_with_sparse_grid(double (*function)(std::vector<double> x, Args... rest),int d, int l, std::vector<std::vector<double> > nodes,std::vector<std::vector<double> > weights,bool write_in_file, bool use_trap_rule,  Args... rest){
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
			printf("%i ki \n", k[i]);
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

		for(int i=0; i<d; i++){
		trap_rule_s(&weights[i], k[i]);
		trap_rulen(&nodes[i], k[i]);
		}
		std::cout<<nodes[1][2]<<" nodes "<<std::endl;
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
				std::cout<<point[h]<<" point "<<std::endl;
			}
			final_value = final_value + product_w*(function_to_integrate(point, rest...));
			product_w = 1;
	//		std::cout<<product<<std::endl;
	//		std::cout<<final_value<<std::endl;

		}

	}

	return final_value;

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




/*int main(){


	int d=2;
	int l=7;
	int maxlevel = pow(2,l)-1;



	double** nodes = new double*[d];
	double** weights = new double*[d];
	for(int i=0; i<d; i++) {
		weights[i] =  new double[maxlevel];
		nodes[i] = new double[maxlevel];
	}

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

	double final_value = integrate_with_sparse_grid(testfunction,
			d, l, nodesv, weightsv);

	printf("%f final sum\n", final_value);


}
*/



//////////////////////////////////////////////////////////
//////////////////////////Task_2//////////////////////////
//////////////////////////////////////////////////////////
double discrete_geometric_average_exact(double s0, double r, double T, int M, double K, double sigma)
{
	double delta_t = T/(double)M;
	double T1 = T-((M*(M-1.)*(4.*M+1.))/(6.*M*M))*delta_t;
	double T2 = T-((M-1.)/2.)*delta_t;
	double A = exp(-r*(T-T2)-(sigma*sigma*(T2-T1))/2);
	double d = (log(s0/K)+(r-0.5*sigma*sigma)*T2)/(sigma*sqrt(T1));

	return s0*A*normal_cdf(d+sigma*sqrt(T1))-K*exp(-r*T)*normal_cdf(d);
}

double continuous_geometric_average_exact(double s0, double r, double T, double K, double sigma)
{
	double d = (log(s0/K)+0.5*(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T/3.));

	return s0*exp(-0.5*(r+(sigma*sigma)/6.)*T)*normal_cdf(d+sigma*sqrt(T/3.))-K*exp(-r*T)*normal_cdf(d);
}

double discrete_geometric_average_simulation(gsl_rng* rng, double s0, double r, double T, int M, double K, double sigma, int N)
{
	double delta_t = T/(double)M;
	double result = 0.0;
	std::vector<double>* w;
	std::vector<double>* s;

	for(int i=0; i<N; i++)
	{
		/* Simulation of brownian motion */
		w = wiener_process(rng, T, delta_t);
		s = brownian_motion(rng, T, delta_t, w, s0, r, sigma);

		double product = 1.0;
		for(int j=0; j<M; j++)
		{
			product *= (*s)[j];
		}
		product = pow(product, 1./(double)M)-K;
		if(product > 0)
			result += product;
	}
	result = result/(double)N;
	return result;
}

//////////////////////////////////////////////////////////
//////////////////////////Task_8//////////////////////////
//////////////////////////////////////////////////////////

/*
double payoff_discrete_arithmetic_average(gsl_rng* rng, double s0, double r, double T, int M, double K, double sigma)
{
	double sum = 0.0;
	double delta_t = T/(double)M;
	std::vector<double>* w;
	std::vector<double>* s;

	// Simulation of brownian motion
	w = wiener_process(rng, T, delta_t);
	s = brownian_motion(rng, T, delta_t, w, s0, r, sigma);

	for(int i=0; i<M; i++)
	{
		sum += (*s)[i];
	}
	sum = sum/(double)M;

	if(sum-K<0.0)
		return 0.0;
	else
		return sum-K;
}
*/


double asian_option_call_integrand(std::vector<double> x,double S0,double K, double sigma, double mu,int M, double T, bool bb_not_rw){
  // computing num of levels;
  int max_level = log2(M);
  double result = S0;
  double delta_t = T/M;

  std::vector<double>* w = new std::vector<double>(M);
  std::vector<double>::iterator it = w->begin();
  w->clear();

  if (bb_not_rw==false) {
    (*w)[0] = 0.0;
    for (int i = 1; i < M; i++) {
      result *= S0*exp((mu-0.5*sigma*sigma)*i*delta_t+sigma*((*w)[i-1]+delta_t*x[i]));
      (*w)[i] = (*w)[i-1]+delta_t*x[i];
    }
  }
  else if(bb_not_rw==true){
    for (size_t i = 1; i <= max_level; i++) {
      for (size_t j = 0; j < pow(2,i); j++) {
        result *= S0*exp((mu-0.5*sigma*sigma)*i*delta_t+sigma*(0.5*((*w)[i]+(*w)[i+1])+delta_t*x[i]));
        w->insert(it+i+1, 0.5*((*w)[i]+(*w)[i+1])+delta_t*x[i]);
      }
    }
  };

  result = pow(result,1./M)-K;
  if(result>0.0){
    return result;
  }
  else{
    return 0.0;
  }


  }


  //////////////////////////////////////////////////////////
  //////////////////////////Task_8//////////////////////////
  //////////////////////////////////////////////////////////
  /* Nl has all the different N_l_1,....,N_l_d */
  double function_to_integrate(std::vector<double> x)
  {
  return x[0]+x[1]+x[2];
  }


  //////////////////////////////////////////////////////////
  //////////////////////////Task_9//////////////////////////
  //////////////////////////////////////////////////////////
  template<typename... Args>
  void write_quadrature_points_to_file(std::ofstream& myfile, int iteration, std::vector<std::vector<double>> nodes_temp, int d, std::vector<int> Nl, std::vector<int> ids)
  {
  	if(iteration==d)
  	{
  		for(int i=0; i<d; i++)
  		{
  			myfile << nodes_temp[i][ids[i]] << "	";
  			//std::cout << nodes_temp[i][ids[i]] << "	";
  		}
  		myfile << std::endl;
  		//std::cout << std::endl;
  	}
  	else
  	{
  		for(int k=0; k<Nl[iteration]; k++)
  		{
  			ids[iteration] = k;
  			write_quadrature_points_to_file(myfile, iteration+1, nodes_temp, d, Nl, ids);
  		}
  	}
  }


  //////////////////////////////////////////////////////////
  //////////////////////////Task_6//////////////////////////
  //////////////////////////////////////////////////////////
  std::vector<double> van_der_corput_sequence(int p, int n, double epsilon)
  {
  	std::vector<double> x;
  	double x_i_1 = 0;
  	double z;
  	double v;

  	for(int i=0; i<n; i++)
  	{
  		z = 1-x_i_1;
  		v=1./p;
  		while(z<v+epsilon)
  		{
  			v=v/p;
  		}
  		x_i_1 = x_i_1+(p+1.)*v-1.;
  		x.push_back(x_i_1);
  	}

  	return x;
  }

  bool is_prime(int number)
  {
  	for(int i=2; i<=sqrt(number); i++)
  	{
  		if(number%i == 0)
  			return false;
  	}
  	return true;
  }

  std::vector<int> first_prime_numbers(int n)
  {
  	std::vector<int> prime_numbers;
  	int i = 0;
  	int number = 2;

  	while(i<n)
  	{
  		if(is_prime(number))
  		{
  			prime_numbers.push_back(number);
  			i++;
  		}
  		number++;
  	}

  	return prime_numbers;
  }

  std::vector<std::vector<double>> d_dimensional_halton_sequence(int d, int n)
  {
  	std::vector<std::vector<double>> points;
  	std::vector<int> prime_numbers = first_prime_numbers(d);

  	for(int i=0; i<n; i++)
  	{
  		std::vector<double> single_point;
  		points.push_back(single_point);
  	}

  	for(int j=0; j<d; j++)
  	{
  		std::vector<double> van_der_corput_sequence_j = van_der_corput_sequence(prime_numbers[j], n, pow(10.,-12.));
  		for(int i=0; i<n; i++)
  		{
  			points[i].push_back(van_der_corput_sequence_j[i]);
  		}
  	}

  	return points;
  }

  std::vector<double>* brownian_bridge_level(std::vector<double>* prev_level, double T, double level, gsl_rng* r){

  	//definition of iterator
  	std::vector<double>::iterator it = prev_level->begin();

  	double new_point = 0;
  	for (int i = 0; i < pow(2,level); i+=2) {
  		//generating points of new level from thr previous one and inserting them inbetween
  		new_point = (1./2.) * ( (*prev_level)[i] + (*prev_level)[i+1]) + sqrt(T/pow(2,level))*gsl_ran_ugaussian(r);
  		prev_level->insert(it+i+1,new_point);

  	}
  	// returns pointer to the vector with points of next level
  	return prev_level;
  }

  std::vector<double>* brownian_bridge(gsl_rng* r, double T, int M)
  {
  		// computing number of levels
  		int max_level = log2(M);

  		// start settings
  		std::vector<double> *w = new std::vector<double>;
  		w->clear();
  		w->push_back(0.0);
  		w->push_back(gsl_ran_ugaussian(r));

  		// recursion similar call for computation of new levels
  		for(int i=1; i<=max_level; i++)
  		{
  			brownian_bridge_level(w,T,i,r);
  		}
  		return w;
  }


  void monte_carlo_multivariate(std::vector<std::vector<double>>* nodes, std::vector<double>* weights, int l, int d, gsl_rng* r)
  {
  	int Nl = pow(2,l)-1;
  	for(int i=1; i<=Nl; i++)
  	{
  		for(int j=1;j<=d;j++){
  		(*nodes)[i-1].push_back(random_number_01_GSL(r));
  		}
  		(*weights)[i-1] = (1./(Nl));
  	}
  }

 template<typename... Args>
  double integrate_by_point_evaluation_multivariate(double (*function)(std::vector<double>* x, Args... rest),int n, int d, std::vector<std::vector<double>>* nodes, std::vector<double>* weights, Args... rest)
  {
  	double result = 0.0;
  	for(int i=0; i<n; i++)
  	{
  		result += (*weights)[i]*function(&(*nodes)[i], rest...);
  	}

  	return result;
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