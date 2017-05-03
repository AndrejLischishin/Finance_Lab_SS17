#include<cstdlib>
#include<cstdio>
#include<gsl/gsl_rng.h>
#include<math.h>

#define A0 0.398942270991
#define A1 0.020133760596
#define A2 0.002946756074
#define B1 0.217134277847
#define B2 0.018576112465
#define B3 0.000643163695
#define C0 1.398247031184
#define C1 -0.360040248231
#define C2 0.022719786588
#define D0 1.460954518699
#define D1 -0.305459640162
#define D2 0.038611796258
#define D3 -0.003787400686
#define E0 2.50662823884
#define E1 -18.61500062529
#define E2 41.39119773534
#define E3 -25.44106049637
#define F0 -8.47351093090
#define F1 23.08336743743
#define F2 -21.06224101826
#define F3 3.13082909833
#define G0 0.3374754822726147
#define G1 0.9761690190917186
#define G2 0.1607979714918209
#define G3 0.0276438810333863
#define G4 0.0038405729373609
#define G5 0.0003951896511919
#define G6 0.0000321767881768
#define G7 0.0000002888167364
#define G8 0.0000003960315187

// Task 1

double randomNumber01()
{
	return (double) rand() / RAND_MAX; // *
}

double randomNumber01GSL(gsl_rng* r)
{
	return gsl_rng_uniform(r);
}

// Task 2

double randomNumberabGSL(gsl_rng* r, double a, double b)
{	
	return a + gsl_rng_uniform(r)*(b-a);	
}

double rej_sam_alg(double a, double b, double max, gsl_rng* r){
	double x = randomNumberabGSL(r, a, b);
	double y = randomNumberabGSL(r, 0, max);
	if(y<=1./sqrt(2*M_PI)*exp(-x*x*0.5)){
		return x;
	}
	else rej_sam_alg(a, b, max, r);
	
	return 0;
}

// Task 4

double normal_CDF(double x){
	double x2;
	if(x<0.0) return 1.0-normal_CDF(-x);
	if(x<=1.87) {
	x2=x*x;
	return 0.5+x*(A0+(A1+A2*x2)*x2)/(1.0+(B1+(B2+B3*x2)*x2)*x2);
	}
	else if(x<6.0){
	return 1.0 - pow((C0+(C1+C2*x)*x)/(D0+(D1+(D2+D3*x)*x)*x), 16.0);
	}
	return 1.0;
}

double normal_inverse_CDF(double x){
    double p = x-.5;
    double r;
    if(0.42>fabs(p)){
        r = p*p;
        return p*(((E3*r+E2)*r+E1)*r+E0)/((((F3*r+F2)*r+F1)*r+F0)*r+1.0);
    }
    else{
        if(p<0) r = x;
        else r = 1 - x;
        
        r = log(-log(r));
        r = G0+r*(G1+r*(G2+r*(G3+r*(G4+r*(G5+r*(G6+r*(G7+r*G8)))))));

        if(p<0) return -r;
        else    return r;
    }
}

double normal_inverse_CDF_ex_var(double x, double ex, double sigma){
    double p = x-.5;
    double r;
    if(0.42>fabs(p)){
        r = p*p;
        return p*(((E3*r+E2)*r+E1)*r+E0)/((((F3*r+F2)*r+F1)*r+F0)*r+1.0);
    }
    else{
        if(p<0) r = x;
        else r = 1 - x;
        
        r = log(-log(r));
        r = G0+r*(G1+r*(G2+r*(G3+r*(G4+r*(G5+r*(G6+r*(G7+r*G8)))))));

        if(p<0)  return ex-r*sigma;
        else     return ex+r*sigma;
    }
}

// Task 6

void box_muller(double* z_1, double* z_2, int n, gsl_rng* r){
	for(int i=0; i<n; i++){
		double u_1 = randomNumber01GSL(r);
		double u_2 = randomNumber01GSL(r);
		z_1[i]=sqrt(-2*log(u_1))*cos(2*M_PI*u_2);
		z_2[i]=sqrt(-2*log(u_1))*sin(2*M_PI*u_2);
	}
}

void box_muller_ex_var(double* z_1, double* z_2, int n, gsl_rng* r, double x, double sigma){
	for(int i=0; i<n; i++){
		double u_1 = randomNumber01GSL(r);
		double u_2 = randomNumber01GSL(r);
		z_1[i]=x + sqrt(-2*log(u_1))*cos(2*M_PI*u_2)*sigma;
		z_2[i]=x + sqrt(-2*log(u_1))*sin(2*M_PI*u_2)*sigma;
	}
}

// Task 8

double var_est(double* values, int n){
	double alpha = values[0];
	double beta = 0;
	double sigma;
	
	for(int i=0; i<n; i++){
		double gamma = values[i] - alpha;
		alpha = alpha + gamma/(i+1);
		beta = beta+gamma*gamma*i/(i+1);
		sigma = sqrt(beta/i);
	}
	return sigma;	
}

 // Task10
 
void wiener_process(double* w, double ex, double sigma, double delta_t, double T, gsl_rng* r){
	
	w[0]=0;
	double M = (double) T/delta_t;
	for(int i=0; i<=M; i++){
	
	double x = randomNumber01GSL(r);
	x = normal_inverse_CDF_ex_var(x, ex, sigma);
	
	w[i+1] = w[i] + sqrt(delta_t)*x;
	}	
}