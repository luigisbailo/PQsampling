#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>

#include "../parameters.hpp"
#include "../draw.hpp"

int main () {


	double b=1;
	double D=1;
	double tau = 1;
 	double t1 = 0.01;
 	double t2 = 0.1;
 	double t3 = 0.99;
 	double r;

	PQ00der (0.1,0.1,tau,b,D,Sder(tau,b,D));

	for ( double r=b/100; r<b; r+=b/100){


		std::cout << r << "\t"
			 << Pder(r,t1,b,D,Sfunct(t1,b,D)) << "\t"<< Pder(r,t2,b,D,Sfunct(t2,b,D)) << "\t"<< Pder(r,t3,b,D,Sfunct(t3,b,D)) << "\t"
			 << PQ00der (r,t1,tau,b,D,Sder(tau,b,D)) << "\t" << PQ00der (r,t2,tau,b,D,Sder(tau,b,D)) << "\t" <<PQ00der (r,t3,tau,b,D,Sder(tau,b,D)) << std::endl;

	}




}
