//
//  exotic_options.cpp
//

#include "exotic_options.hpp"

double payoff_discrete_down_out_call(std::vector<double>* x, double s0, double r, double T, int M, double K, double sigma, double B)
{
	std::vector<double> z;
	std::vector<double> S;
    std::vector<double> w(M+1);
	double delta_t = (double)T/M;
    w[0] = 0;
    
    for(int i=0; i<M; i++)
	{
		z.push_back(normal_inverse_cdf((*x)[i]));
        w[i+1] = w[i]+sqrt(delta_t)*z[i];
		S.push_back(s0*exp((r-(sigma*sigma)/2.0)*delta_t*(i+1.0)+sigma*(w[i+1])));
        if(S[i]<=B)
			return 0.0;
	}

	double temp = S[M-1]-K;

	if(temp>0.0)
		return temp/exp(r*T);
	else
		return 0.0;
}

double payoff_discrete_lookback(std::vector<double>* x, double s0, double r, double T, int M, double K, double sigma)
{
	std::vector<double> z;
	std::vector<double> S;
    std::vector<double> w(M+1);
	double delta_t = (double)T/M;
    w[0] = 0;

	for(int i=0; i<M; i++)
	{
		z.push_back(normal_inverse_cdf((*x)[i]));
        w[i+1] = w[i]+sqrt(delta_t)*z[i];
		S.push_back(s0*exp((r-(sigma*sigma)/2.0)*delta_t*(i+1.0)+sigma*w[i+1]));
        
	}

	double max = S[0];

	for(int i=1; i<M; i++)
	{
		if(S[i]>max)
		{
			max = S[i];
		}
	}

	double temp = max-K;

	if(temp > 0)
		return temp/exp(r*T);
	else
		return 0.0;
}


//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////barrier_integrand_fixed/\/\bb/\or not////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double barrier_integrand(std::vector<double>* x,double S0,double K, double sigma, double mu,int M, double T, double B, bool use_bb){
    // computing num of levels;
    int max_level = log2(M);
    double result = 0;//S0;
    
    double delta_t = 0;
    double helper = 0;
    bool tracker = true;
    
    
    std::vector<double> w(2*(M+1),0);
    
    if (use_bb==false) {
        delta_t = T/M;
        w[0] = 0.0;
        for (int i = 1; i <= M; i++) {
            result = S0*exp((mu-0.5*sigma*sigma)*i*delta_t+sigma*(w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1])));
            w[i] = w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1]);
            if (result<=B) {
                tracker = false;
                break;
            }
        }
        w.resize(M);
    }
    else if(use_bb==true){
        w[0] = 0.0;
        w[1]=sqrt(T)*normal_inverse_cdf((*x)[M-1]);
        int count = 0;
        
        for (int i = 1; i <= max_level; i++) {
            delta_t = T/(pow(2,(double)i));
            for (int j = 0; j < pow(2,i); j+=2) {
                helper = 0.5*(w[j]+w[j+1])+sqrt(delta_t/2.)*normal_inverse_cdf((*x)[count]);
                result = S0*exp((mu-0.5*sigma*sigma)*(j+1)*delta_t+sigma*(helper));
                w.emplace(w.begin()+j+1, helper);
                count++;
                if (result<=B) {
                    tracker = false;
                    break;
                }
                
            }
            if (tracker == false) {
                break;
            }
        }
        w.resize(M);
        //result = result*trsult(T)
        if (tracker == true) {
            result = S0*exp((mu-0.5*sigma*sigma)*T+sigma*(sqrt(T)*normal_inverse_cdf((*x)[M-1])));
            if (result<=B) {
                tracker = false;
            }
        }
        
    };
    
    
    if (tracker==true) {
        result = exp(-mu*T)*(result-K);
        if(result>0.0){
            return result;
        }
        else return 0.0;
    }
    else{
        return 0.0;
    };
    
}


//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////look_back_integrand_fixed/\/\bb/\or not////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double lookback_call_integrand_fixed(std::vector<double>* x,double S0,double K, double sigma, double mu,int M, double T, bool use_bb){
    // computing num of levels;
    int max_level = log2(M);
    double result = 0;//S0;
    double max_s = 0;
    double delta_t = 0;
    double helper = 0;
    
    
    std::vector<double> w(2*(M+1),0);
    
    if (use_bb==false) {
        delta_t = T/M;
        w[0] = 0.0;
        for (int i = 1; i <= M; i++) {
            result = S0*exp((mu-0.5*sigma*sigma)*i*delta_t+sigma*(w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1])));
            w[i] = w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1]);
            if (result>max_s) {
                max_s = result;
            }
        }
        w.resize(M);
    }
    else if(use_bb==true){
        w[0] = 0.0;
        w[1]=sqrt(T)*normal_inverse_cdf((*x)[M-1]);
        int count = 0;
        
        for (int i = 1; i <= max_level; i++) {
            delta_t = T/(pow(2,(double)i));
            for (int j = 0; j < pow(2,i); j+=2) {
                helper = 0.5*(w[j]+w[j+1])+sqrt(delta_t/2.)*normal_inverse_cdf((*x)[count]);
                result = S0*exp((mu-0.5*sigma*sigma)*(j+1)*delta_t+sigma*(helper));
                w.emplace(w.begin()+j+1, helper);
                count++;
                if (result>max_s) {
                    max_s = result;
                }
                
            }
        }
        w.resize(M);
        //result = result*trsult(T)
        result = S0*exp((mu-0.5*sigma*sigma)*T+sigma*(sqrt(T)*normal_inverse_cdf((*x)[M-1])));
        if (result>max_s) {
            max_s = result;
        }
    };
    
    
    result = exp(-mu*T)*(max_s-K);
    if(result>0.0){
        return result;
    }
    else{
        return 0.0;
    };
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////asian_option_call_integrand_arithmetic/\/\/\for geometric see multivariateintegration.cpp///
///////////////////////////////////////////////////////////////////////////////////////////////////


double asian_option_call_integrand_arithmetic(std::vector<double>* x,double S0,double K, double sigma, double mu,int M, double T, bool use_bb){
    // computing num of levels;
    int max_level = log2(M);
    double result = 0;//S0;
    double delta_t = 0;
    double helper = 0;
    
    
    std::vector<double> w(2*(M+1),0);
    
    if (use_bb==false) {
        delta_t = T/M;
        w[0] = 0.0;
        for (int i = 1; i <= M; i++) {
            result += S0*exp((mu-0.5*sigma*sigma)*i*delta_t+sigma*(w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1])));
            w[i] = w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1]);
        }
        w.resize(M);
    }
    else if(use_bb==true){
        w[0] = 0.0;
        w[1]=sqrt(T)*normal_inverse_cdf((*x)[M-1]);
        int count = 0;
        
        for (int i = 1; i <= max_level; i++) {
            delta_t = T/(pow(2,(double)i));
            for (int j = 0; j < pow(2,i); j+=2) {
                helper = 0.5*(w[j]+w[j+1])+sqrt(delta_t/2.)*normal_inverse_cdf((*x)[count]);
                result += S0*exp((mu-0.5*sigma*sigma)*(j+1)*delta_t+sigma*(helper));
                w.emplace(w.begin()+j+1, helper);
                count++;
                
                
            }
        }
        w.resize(M);
        //result = result*trsult(T)
        result = result + S0*exp((mu-0.5*sigma*sigma)*T+sigma*(sqrt(T)*normal_inverse_cdf((*x)[M-1])));
        
    };
    
    
    result = exp(-mu*T)*((result/(double)M)-K);
    if(result>0.0){
        return result;
    }
    else{
        return 0.0;
    };
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////asian_option_call_integrand_arithmetic/\/\/\for geometric see multivariateintegration.cpp///
///////////////////////////////////////////////////////////////////////////////////////////////////


double asian_option_call_integrand_arithmetic_control_variates(std::vector<double>* x,double S0,double K, double sigma, double mu,int M, double T, bool use_bb){
    // computing num of levels;
    int max_level = log2(M);
    double result_arithmetic = 0;//S0;
    long double result_geometric = 1;
    double result;
    double delta_t = 0;
    double helper = 0;
    
    
    
    std::vector<double> w(2*(M+1),0);
    
    if (use_bb==false) {
        delta_t = T/M;
        w[0] = 0.0;
        for (int i = 1; i <= M; i++) {
            result_arithmetic += S0*exp((mu-0.5*sigma*sigma)*i*delta_t+sigma*(w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1])));
            result_geometric *= S0*exp((mu-0.5*sigma*sigma)*i*delta_t+sigma*(w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1])));
            w[i] = w[i-1]+sqrt(delta_t)*normal_inverse_cdf((*x)[i-1]);
        }
        w.resize(M);
    }
    else if(use_bb==true){
        w[0] = 0.0;
        w[1]=sqrt(T)*normal_inverse_cdf((*x)[M-1]);
        int count = 0;
        
        for (int i = 1; i <= max_level; i++) {
            delta_t = T/(pow(2,(double)i));
            for (int j = 0; j < pow(2,i); j+=2) {
                helper = 0.5*(w[j]+w[j+1])+sqrt(delta_t/2.)*normal_inverse_cdf((*x)[count]);
                result_arithmetic += S0*exp((mu-0.5*sigma*sigma)*(j+1)*delta_t+sigma*(helper));
                result_geometric *= S0*exp((mu-0.5*sigma*sigma)*(j+1)*delta_t+sigma*(helper));
                w.emplace(w.begin()+j+1, helper);
                count++;
                
                
            }
        }
        w.resize(M);
        //result = result*trsult(T)
        result_arithmetic += S0*exp((mu-0.5*sigma*sigma)*T+sigma*(sqrt(T)*normal_inverse_cdf((*x)[M-1])));
        result_geometric *=  S0*exp((mu-0.5*sigma*sigma)*T+sigma*(sqrt(T)*normal_inverse_cdf((*x)[M-1])));
        
    };
    
    
    result_arithmetic = exp(-mu*T)*((result_arithmetic/(double)M)-K);
    result_geometric = exp(-mu*T)*(pow(result_geometric, 1./(double)M)-K);
    if (result_arithmetic < 0.0) {
        result_arithmetic = 0.0;
    };
    
    if (result_geometric < 0.0) {
        result_geometric = 0.0;
    };
    
    result = result_arithmetic - result_geometric;
    return result;
}

double black_scholes_down_out_call(double s0, double K, double T, double sigma, double r, double B)
{
	double Z;
	if(B==0)
		Z = 0;
	else
		Z = pow((B/s0),(((2.0*r)/(sigma*sigma))-1.0));

	double Bbar;

	if(B>K)
		Bbar = B;
	else
		Bbar = K;

	return V_bs(s0, Bbar, r, sigma, T)-Z*V_bs(B*B/s0, Bbar, r, sigma, T)+(Bbar-K)*exp(-r*T)*(normal_cdf(d_S_K(s0, Bbar, r, sigma, T))-Z*normal_cdf(d_S_K(B*B/s0, Bbar, r, sigma, T)));
}

double d_S_K(double S, double K, double r, double sigma, double T)
{
	return (log(S/K)+(r-sigma*sigma/2.0)*T)/(sigma*sqrt(T));
}

double V_bs(double S, double K, double r, double sigma, double T)
{
	double d = d_S_K(S, K, r, sigma, T);
	return S*normal_cdf(d+sigma*sqrt(T))-K*exp(-r*T)*normal_cdf(d);

}
