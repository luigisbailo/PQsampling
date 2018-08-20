//TO COMPILE: g++ -std=c++11 main.cpp -o main -lgsl -lgslcblas -lm


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <sstream>


#include "../parameters.hpp"
#include "../draw.hpp"
#include "../tools.hpp"
#include "../init.hpp"
#include "../step.hpp"
#include "../shell.hpp"
#include "../print.hpp"
#include "../burst.hpp"
#include "../bruteForce.hpp"
#include "../checks.hpp"
#include "../run_BM.hpp"


int main (int argc, char *argv[]) {


	
	double D_A = 0.01;
	double D_B = 0.01;
	double R_A = 2.5;
	double R_B = 2.5;
	
	double tau_bm = 0.1;
	const int N = 10; 
	const int N_A = 5;
	const int N_B = 5;

	double alpha= 9;
	double L = 20;
	// int Nsamples = 100;

	int nProj;
	double Tsim;
	int Nsamples;

	std::stringstream convert_nProj (argv[1]);
  	std::stringstream convert_Tsim (argv[2]);
  	std::stringstream convert_Nsamples (argv[3]);

	if (!(convert_nProj >> nProj ))
	    exit (EXIT_FAILURE);  
	if (!(convert_Tsim >> Tsim ))
	    exit (EXIT_FAILURE);  
	if (!(convert_Nsamples >> Nsamples ))
	    exit (EXIT_FAILURE);  

	int stat[3];
	// double Diff_BM [N][Nsamples][nProj];
	// double diffStat[nProj][10];
	double Diff_BM [nProj];
	double diffStat[nProj];


	for ( int t=0; t<nProj; t++)
		Diff_BM[t] = 0;


	for ( int count = 0; count < Nsamples; count++){


			for ( int t=0; t<nProj; t++)
				diffStat[t] = 0;

		run_BM ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, Tsim, nProj, L, diffStat );


		for ( int t=0; t<nProj; t++){

			Diff_BM [t] += diffStat[t];

		}



	}	


	for ( int t=0; t<nProj; t++){

		std::cout << (t+1)*Tsim/nProj << "\t" << Diff_BM[t]/(Nsamples*N) <<  std::endl ;

	}


}
