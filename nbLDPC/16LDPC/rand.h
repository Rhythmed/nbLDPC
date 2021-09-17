#ifndef _GAUSS_RANDOM_
#define _GAUSS_RANDOM_

#include <math.h>
#include <stdlib.h>


double GaussRand(double dExpect = 0.0, double dVariance = 1.0);
double GaussRand(double dExpect, double dVariance){
	static double V1, V2, S = 0.0;
	static int phase = 0;
	double X = 0.0;

	if (phase == 0){
		do{
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);
		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;
	return (X*dVariance + dExpect);
}


#endif