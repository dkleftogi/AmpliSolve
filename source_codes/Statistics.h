#ifndef STATISTICS_H_INCLUDED
#define STATISTICS_H_INCLUDED

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "DataStructures.h"

/////////////////////////////////////////////////////////////////////////////////////
// Normal CDF
/////////////////////////////////////////////////////////////////////////////////////

double erf(double x)
{
     double y = 1.0 / ( 1.0 + 0.3275911 * x);
     return 1 - (((((
            + 1.061405429  * y
            - 1.453152027) * y
            + 1.421413741) * y
            - 0.284496736) * y
            + 0.254829592) * y)
            * exp(-x * x);
}

double pdf(double x, double mu, double sigma)
{

      const double pi = 3.14159265;
      return exp( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
}

double cdf(double x, double mu, double sigma)
{
	return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2.))));
}

/////////////////////////////////////////////////////////////////////////////////////
// Fisher test
/////////////////////////////////////////////////////////////////////////////////////

void initLogFacs(double* logFacs , int n)
{
    int i;
    logFacs[0] = 0;
    for(i=1; i < n+1; ++i)
    {
        logFacs[i] = logFacs[i-1] + log((double)i);
    }
}

double logHypergeometricProb(double *logFacs, int a, int b, int c, int d)
{
    return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

double fexactt(int a, int b, int c, int d)
{
    int n = a + b + c + d;

    double *logFacs = (double*)malloc(sizeof(double)*(n+1));
    initLogFacs(logFacs,n);

    double logpCutoff = logHypergeometricProb(logFacs,a,b,c,d);
    double pFraction = 0;

    int x;

    for(x=0;x<=n;++x)
    {
        if (a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0) {
            double l = logHypergeometricProb(logFacs,x,a+b-x,a+c-x,d-a+x);
            if (l <= logpCutoff)
                pFraction += exp(l-logpCutoff);
        }
    }

    double logpValue = logpCutoff + log(pFraction);

    return(exp(logpValue));
}

#endif // STATISTICS_H_INCLUDED
