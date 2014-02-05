#ifndef RANDOM_H
#define RANDOM_H

#include <vector>

double ran2(long *);
double ran_unif(double , double);
// Bernoulli(0,1), p is the probability for success.
inline bool ran_ber(double p) {
	return ran_unif(0, 1) <= p;
}
// Integrated uniform random number(n,m)
inline long ran_iunif(long min, long max) {
	return (long)ran_unif((double)min, (double)(max+1));
}
// Return specified two numbers a or b
// a indicates the success one with the probability p.
inline long ran_nber(long a, long b, double p) {
	return ran_ber(p)? a: b;
}



double ran_norm(double , double );
double ran_exp(double lambda);
double ran_gamma(double , double );
long ran_num(const std::vector<double> &p, const std::vector<long> &num);
unsigned long long ran_num(const std::vector<double> &p, 
	const std::vector<unsigned long long> &num);
unsigned long long ran_num_log(const std::vector<double> &p, const std::vector<unsigned long long> &num);
double gammaln(double xx);
double gammaf(double xx);
double betaf(double x, double y);
double betad(double x, double a, double b);

double f_fmin(double x, double y);
double f_fmax(double x, double y);
double ran_beta(double aa, double bb);













#endif