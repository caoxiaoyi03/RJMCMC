#include"Random.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <float.h>
#include <vector>
#include <map>
#include <algorithm>
using namespace std;


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define E 2.71828183

double ran2(long *idum)
{
	long j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



// Uniform(min,max)
double ran_unif(double min, double max) {
	static long idum = -1;
	double rtemp;
	if(min >= max) return 0;
	if(idum < 0) {
		// not initialized
		srand((unsigned)time(NULL));
		idum *= rand();
		while(!idum)
			idum = -rand();
	}
	do {
		rtemp = ran2(&idum);
		rtemp = rtemp * (max - min) * 2 + min * 1.5 - max * 0.5;
	} while (rtemp <= min || rtemp >=max);

	return rtemp;

}

////////////////////////////////////
//Normal(u,sigma)double ran_norm(double , double );
double ran_norm(double mu, double sigma)
{

	static long idum = -1;
	if(idum < 0) {
		// not initialized
		srand((unsigned)time(NULL));
		idum *= rand();
		while(!idum)
			idum = -rand();
	}
	static long iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran2(&idum)-1.0;
			v2=2.0*ran2(&idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return (v2*fac)*sigma+mu;
	} else {
		iset=0;
		return (gset)*sigma+mu;
	}
}


// Exponential random number
double ran_exp(double lambda)
{
 
	static long idum = -1;
	if(idum < 0) {
		// not initialized
		srand((unsigned)time(NULL));
		idum *= rand();
		while(!idum)
			idum = -rand();
	}

	double dum;

	do
		dum=ran2(&idum);
	while (dum == 0.0);
	if(lambda>0)
	{
	 return -log(dum)/lambda;
	} else
	{
	 cout<<"Error! Lambda is not positive!"<<endl;
	 return 0.0;
	}
}




// Gamma: Debasis Kundu (2007): A Convenient Way of Generating Gamma
//Random Variables Using Generalized
//Exponential Distribution, algorithm 3

#define repeat for(;;)
double ran_gamma(double shape, double scale)
{
 // First generate gamma(shape, scale=1)
 // For 0<shape<1
	   double d = 1.0334 - 0.0766*pow(E,2.2942*shape);
       double a = pow(2,shape)*pow(1-pow(E,-d/2),shape);
       double b = shape*pow(d,shape-1)*pow(E,-d);
       double c = a + b;
       double u,v,x=0;
	   double e, q, r, t, w, ret_val;
       const static double sqrt32 = 5.656854;
    const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

    /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
     * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
     * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
     */
    const static double q1 = 0.04166669;
    const static double q2 = 0.02083148;
    const static double q3 = 0.00801191;
    const static double q4 = 0.00144121;
    const static double q5 = -7.388e-5;
    const static double q6 = 2.4511e-4;
    const static double q7 = 2.424e-4;

    const static double a1 = 0.3333333;
    const static double a2 = -0.250003;
    const static double a3 = 0.2000062;
    const static double a4 = -0.1662921;
    const static double a5 = 0.1423657;
    const static double a6 = -0.1367177;
    const static double a7 = 0.1233795;

    /* State variables [FIXME for threading!] :*/
    static double aa = 0.;
    static double aaa = 0.;
    static double s, s2, dd;    /* no. 1 (step 1) */
    static double q0, bb, si, cc;/* no. 2 (step 4) */

	if((shape>0)&&(shape<1)){

 
       u = ran_unif(0,1);
       if(u<=(a/(a+b)))
	       x = -2.0*log(1-pow(c*u,1/shape)/2);
       else
	       x = -1.0*log(c*(1-u)/(shape*pow(d,shape-1)));
       v = ran_unif(0,1);
       while(((x>d)||(v>(pow(x,shape-1)*pow(E,-x/2))/(pow(2,shape-1)*pow(1-pow(E,-x/2),shape-1))))&&((x<=d)||(v>pow(d/x,1-shape))))
	   {
        u = ran_unif(0,1);
        if(u<=(a/(a+b)))
	       x = -2.0*log(1-pow(c*u,1/shape)/2);
        else
	       x = -1.0*log(c*(1-u)/(shape*pow(d,shape-1)));
        v = ran_unif(0,1);
		}
	   return (x*scale);
	} 
	
	
	
	
	
	else if(shape>=1)// --- shape >= 1 : GD algorithm --- 
	{
      	/* Constants : */
    
	
        /* Step 1: Recalculations of s2, s, d if a has changed */
    if (shape != aa) {
	aa = shape;
	s2 = shape - 0.5;
	s = sqrt(s2);
	dd = sqrt32 - s * 12.0;
    }
    /* Step 2: t = standard normal deviate,
               x = (s,1/2) -normal deviate. */

    /* immediate acceptance (i) */
    t = ran_norm(0.0,1.0);
    x = s + 0.5 * t;
    ret_val = x * x;
    if (t >= 0.0)
	return scale * ret_val;

    /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
    u = ran_unif(0.0,1.0);
    if (dd * u <= t * t * t)
	return scale * ret_val;

    /* Step 4: recalculations of q0, b, si, c if necessary */

    if (shape != aaa) {
	aaa = shape;
	r = 1.0 / shape;
	q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
	       + q2) * r + q1) * r;

	/* Approximation depending on size of parameter a */
	/* The constants in the expressions for b, si and c */
	/* were established by numerical experiments */

	if (shape <= 3.686) {
	    bb = 0.463 + s + 0.178 * s2;
	    si = 1.235;
	    cc = 0.195 / s - 0.079 + 0.16 * s;
	} else if (shape <= 13.022) {
	    bb = 1.654 + 0.0076 * s2;
	    si = 1.68 / s + 0.275;
	    cc = 0.062 / s + 0.024;
	} else {
	    bb = 1.77;
	    si = 0.75;
	    cc = 0.1515 / s;
	}
    }
    /* Step 5: no quotient test if x not positive */

    if (x > 0.0) {
	/* Step 6: calculation of v and quotient q */
	v = t / (s + s);
	if (fabs(v) <= 0.25)
	    q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
				      + a3) * v + a2) * v + a1) * v;
	else
	    q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);


	/* Step 7: quotient acceptance (q) */
	if (log(1.0 - u) <= q)
	    return scale * ret_val;
    }

    repeat {
	/* Step 8: e = standard exponential deviate
	 *	u =  0,1 -uniform deviate
	 *	t = (b,si)-double exponential (laplace) sample */
	e = ran_exp(1.0);
	u = ran_unif(0.0,1.0);
	u = u + u - 1.0;
	if (u < 0.0)
	    t = bb - si * e;
	else
	    t = bb + si * e;
	/* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
	if (t >= -0.71874483771719) {
	    /* Step 10:	 calculation of v and quotient q */
	    v = t / (s + s);
	    if (fabs(v) <= 0.25)
		q = q0 + 0.5 * t * t *
		    ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
		      + a2) * v + a1) * v;
	    else
		q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
	    /* Step 11:	 hat acceptance (h) */
	    /* (if q not positive go to step 8) */
	    if (q > 0.0) {
		w = exp(q)-1;
		/* original code had approximation with rel.err < 2e-7 */
		/* if t is rejected sample again at step 8 */
		if (cc* fabs(u) <= w * exp(e - 0.5 * t * t))
		    break;
	    }
	}
    } /* repeat .. until  `t' is accepted */
    x = s + 0.5 * t;
    return scale * x * x;
    }
	
	
	
	
	
	else
	{
	 cout<<"Error! shape parameter is a negative value!";
	 return 0.0;
	}


}
#undef repeat


//Sample a concrete integer with a vector of probabilities
// Note that the length of num should be equal to the p length;
long ran_num(const vector<double> &p, const std::vector<long> &num)
{
 long i,k=0;
 double p_sum = 0;
 long lengthp = p.size();

 map<double, long> pToNum;
 //vector<double> p_sumvector(lengthp);
 double p_sample;

 for(i = 0; i < lengthp; i++)
 {
	 p_sum += p[i];
	 pToNum.insert(map<double, long>::value_type(p_sum, num[i]));
	 //p_sumvector.at(i) = p_sum;
	 
 }
 

 if(p_sum!=0.0){
   return pToNum.lower_bound(ran_unif(0.0, p_sum))->second;
 } else {
  std::cout<<"It is a 0 probability!"<<std::endl;
  return 0;
 }
}


unsigned long long ran_num(const vector<double> &p, 
	const std::vector<unsigned long long> &num)
{
 long i,k=0;
 double p_sum = 0;
 long lengthp = p.size();

 map<double, unsigned long long> pToNum;
 //vector<double> p_sumvector(lengthp);
 double p_sample;

 for(i = 0; i < lengthp; i++)
 {
	 p_sum += p[i];
	 pToNum.insert(map<double, unsigned long long>::value_type(p_sum, num[i]));
	 //p_sumvector.at(i) = p_sum;
	 
 }
 

 if(p_sum!=0.0){
   return pToNum.lower_bound(ran_unif(0.0, p_sum))->second;
 } else {
  std::cout<<"It is a 0 probability!"<<std::endl;
  return 0;
 }
}


//Sample a concrete integer with a vector of probabilities
// Note that the length of num should be equal to the p length;
unsigned long long ran_num_log(const vector<double> &p, const std::vector<unsigned long long> &num)
{
 long i,k=0;
 double p_sum = 0, min = *min_element(p.begin(), p.end());
 long lengthp = p.size();

 map<double, unsigned long long> pToNum;
 //vector<double> p_sumvector(lengthp);
 double p_sample;

 for(i = 0; i < lengthp; i++)
 {
	 p_sum += exp(p[i] - min);
	 pToNum.insert(map<double, unsigned long long>::value_type(p_sum, num[i]));
	 //p_sumvector.at(i) = p_sum;
	 
 }
 

 if(p_sum!=0.0){
   return pToNum.lower_bound(ran_unif(0.0, p_sum))->second;
 } else {
  std::cout<<"It is a 0 probability!"<<std::endl;
  return 0;
 }
}

double gammaln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double gammaf(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return exp(-tmp+log(2.5066282746310005*ser/x));
}

double betaf(double x, double y){
      return gammaf(x)*gammaf(y)/gammaf(x+y);
}


double betad(double x, double a, double b){
   return gammaf(a+b)/gammaf(a)/gammaf(b)*pow(x,a-1)*pow(1-x,b-1);
}

double f_fmin(double x, double y)
{
 return (x>y)?y:x;
}


double f_fmax(double x, double y){

 return (x>y)? x:y;
}



#define M_LN2 0.69314718055994530942
#define expmax	(DBL_MAX_EXP * M_LN2)/* = log(DBL_MAX) */
#define XIFINITE   1999999999

double ran_beta(double aa, double bb)
{
    double a, b, alpha;
    double r, s, t, u1, u2, v, w, y, z;

    int qsame;
    /* FIXME:  Keep Globals (properly) for threading */
    /* Uses these GLOBALS to save time when many rv's are generated : */
    static double beta, gamma, delta, k1, k2;
    static double olda = -1.0;
    static double oldb = -1.0;

    if (aa <= 0. || bb <= 0. || ((aa>XIFINITE) && (bb>XIFINITE))){
	cout<<"ERROR! There are problems in parameter!!!"<<endl;
    return NULL;
	}
    if (aa>XIFINITE)
    	return 1.0;

    if (bb>XIFINITE)
    	return 0.0;

    /* Test if we need new "initializing" */
    qsame = (olda == aa) && (oldb == bb);
    if (!qsame) { olda = aa; oldb = bb; }

    a = f_fmin(aa, bb);
    b = f_fmax(aa, bb); /* a <= b */
    alpha = a + b;

#define v_w_from__u1_bet(AA)			\
	    v = beta * log(u1 / (1.0 - u1));	\
	    if (v <= expmax)			\
		w = AA * exp(v);		\
	    else				\
		w = DBL_MAX                   


    if (a <= 1.0) {	/* --- Algorithm BC --- */

	/* changed notation, now also a <= b (was reversed) */

	if (!qsame) { /* initialize */
	    beta = 1.0 / a;
	    delta = 1.0 + b - a;
	    k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
	    k2 = 0.25 + (0.5 + 0.25 / delta) * a;
	}
	/* FIXME: "do { } while()", but not trivially because of "continue"s:*/
	for(;;) {
	    u1 = ran_unif(0,1);
	    u2 = ran_unif(0,1);
	    if (u1 < 0.5) {
		y = u1 * u2;
		z = u1 * y;
		if (0.25 * u2 + z - y >= k1)
		    continue;
	    } else {
		z = u1 * u1 * u2;
		if (z <= 0.25) {
		    v_w_from__u1_bet(b);
		    break;
		}
		if (z >= k2)
		    continue;
	    }

	    v_w_from__u1_bet(b);

	    if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
		break;
	}
	return (aa == a) ? a / (a + w) : w / (a + w);

    }
    else {		/* Algorithm BB */

	if (!qsame) { /* initialize */
	    beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
	    gamma = a + 1.0 / beta;
	}
	do {
	    u1 = ran_unif(0,1);
	    u2 = ran_unif(0,1);

	    v_w_from__u1_bet(a);

	    z = u1 * u1 * u2;
	    r = gamma * v - 1.3862944;
	    s = a + r - w;
	    if (s + 2.609438 >= 5.0 * z)
		break;
	    t = log(z);
	    if (s > t)
		break;
	}
	while (r + alpha * log(alpha / (b + w)) < t);

	return (aa != a) ? b / (b + w) : w / (b + w);
    }
}

#undef M_LN2
#undef expmax
#undef XIFINITE