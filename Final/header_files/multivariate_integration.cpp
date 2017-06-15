
#include "multivariate_integration.hpp"

void trap_rule_weights(std::vector<double>* weights, int l){
	weights->clear();
	int Nl = pow(2, l)-1;
	int Nk = pow(2, l-1)-1;
	weights->push_back((double)3./(2*(Nl+1)));
	for(int i=1; i<Nl-1; i++){
		 weights->push_back((double)1/(Nl+1));
		 if((i+1)%2==0){
		 	if(i==1 || i==Nl-2){
		 	(*weights)[i] = (*weights)[i]-(double) 3/(2*(Nk+1));
		 	}
		 	else (*weights)[i] = (*weights)[i]-(double) 1/(Nk+1);
		 }
	}
//	(*weights)[0] = (double) 3/(2*(Nl+1));

//	(*weights)[Nl-1] = (double) 3/(2*(Nl+1));
	weights->push_back((double) 3/(2*(Nl+1)));

	if(l==1) (*weights)[0]=1;
	if(l==2) (*weights)[1]=-0.75;
}

void trap_rule_nodes(std::vector<double> *nodes, int l){
	nodes->clear();
	int Nl = pow(2, l)-1;

	for(int i=1; i<=Nl; i++){

		nodes->push_back((double) i/(Nl+1));
	}

}

void clenshaw_curtis_nodes(std::vector<double>* nodes, int l)
{
	nodes->clear();
	unsigned int Nl = pow(2,l)-1;

	for(unsigned int i=1; i<=Nl; i++)
	{
		nodes->push_back(.5*(1.-cos((double)(M_PI*i/(Nl+1.)))));

	}
}


void clenshaw_curtis_weights(std::vector<double>* weights, int l){
	weights->clear();
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
		weights->push_back(( (double)(2./(Nl+1.)) * sin( M_PI * ((double)i/(double)(Nl+1.)) ) * sum ));
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
		}

	}
}


int enumeration(std::vector<int>* k, int d, int l, std::vector<int>* diag){

	int S_k=d;
	int II = 0;
	int counter=0;
	while(1){

	for(int j=0; j<d; j++){
		diag->push_back((*k)[j]-1);
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
			(*k)[j-1]=(*k)[j-1]+1;
			S_k = S_k +1;
			if(S_k > d+l-1){
				if(j==d){ return counter; }
				S_k = S_k-(*k)[j-1]+1;
				(*k)[j-1] = 1;
			}
			else break;
		}

	}

	return counter;
}


void loop(std::vector<int>* vec, std::vector<int>* klevel, int d, std::vector<int>* finalvec){
	finalvec->clear();
	int I;
	int count;
	int number=1;
	int II = 0;
	int last = d-1;

	while(last+1==d){

	for(int j=0; j<d; j++){
//		printf(" %i", vec[j]);
//		finalvec[II+j] = vec[j]-1;
		finalvec->push_back((*vec)[j]-1);
	}
	II = II+d;

//	printf(" end \n");
	number++;
		if((*vec)[last]<(*klevel)[last]){
			(*vec)[last]++;
			}
		else {
		for(int i=d-1; i>=0; i--){
			if((*vec)[i]<(*klevel)[i]) {
			(*vec)[i]++;
			break;
			}

			count = 0;
			for(I=0; I<d; I++){
			if((*vec)[I]==(*klevel)[I]) { count++;
				}
			}
			if(count==d) d=-2;

			if((*vec)[i]==(*klevel)[i]) (*vec)[i]=1;
			}
		}
	}
}

void sparse_grid_nodes(int d, int product, std::vector<int>* allvec, std::vector<std::vector<double> >* nodesv){
	FILE *fp;
	fp = fopen("stuetzstellen", "a");
		for(int i=0; i<product*d; i=i+d){
			for(int j=0; j<d; j++){
				fprintf(fp,"%f ", (*nodesv)[j][(*allvec)[i+j]]);
			}
		fprintf(fp,"\n");
		}
	fclose(fp);
}


double function_task13_sparse(std::vector<double>* x, double gamma, int d){
	double product = 1;
	for(int i=1; i<=d; i++){
		product = product*(1+gamma*exp((*x)[i-1]/2));
	}

	return product;
}


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
		weights->clear();
		nodes->clear();
  	int Nl = pow(2,l)-1;
  	for(int i=1; i<=Nl; i++)
  	{
  		for(int j=1;j<=d;j++){
  		(*nodes)[i-1].push_back(random_number_01_GSL(r));
  		}
			weights->push_back((1./(double)Nl));

  	}
  }

	void quasi_monte_carlo_multivariate(std::vector<std::vector<double>>* nodes, std::vector<double>* weights, int l, int d)
  {
		weights->clear();
		nodes->clear();
  	int Nl = pow(2,l)-1;
  	for(int i=1; i<=Nl; i++)
  	{
  		weights->push_back((1./(double)Nl));
		}
		*nodes = d_dimensional_halton_sequence(d,Nl);
  }



  double function_task13(std::vector<double>* x, double gamma, int d){
  	double product = 1;
  	for(int i=1; i<=d; i++){
  		product = product*(1+gamma*exp((*x)[i-1]/2));
  	}

  	return product;
  }
