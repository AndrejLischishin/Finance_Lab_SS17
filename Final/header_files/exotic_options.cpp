//
//  exotic_options.cpp
//

#include "exotic_options.hpp"

double payoff_discrete_down_out_call(std::vector<double> x, double s0, double r, double T, int M, double K, double sigma, double B)
{
	std::vector<double> z;
	std::vector<double> S;
	double delta_t = (double)T/M;

	for(int i=0; i<M; i++)
	{
		z.push_back(normal_inverse_cdf(x[i]));
		S.push_back(s0*exp((r-(sigma*sigma)/2.0)*delta_t*(i+1.0)+sqrt(delta_t)*z[i]));
		if(S[i]<=B)
			return 0.0;
	}

	double temp = S[M-1]-K;

	if(temp>0.0)
		return temp/exp(r*T);
	else
		return 0.0;
}

double payoff_discrete_lookback(std::vector<double> x, double s0, double r, double T, int M, double K, double sigma)
{
	std::vector<double> z;
	std::vector<double> S;
	double delta_t = (double)T/M;

	for(int i=0; i<M; i++)
	{
		z.push_back(normal_inverse_cdf(x[i]));
		S.push_back(s0*exp((r-(sigma*sigma)/2.0)*delta_t*(i+1.0)+sqrt(delta_t)*z[i]));
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
