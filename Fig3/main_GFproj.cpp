//TO COMPILE: g++ -std=c++11 main.cpp -o main -lgsl -lgslcblas -lm

#define print(x) cout << x << endl;

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

using namespace std::chrono;


#include "../parameters.hpp"
#include "../greensFunct.hpp"
#include "../draw.hpp"
#include "../tools.hpp"
#include "../init.hpp"
#include "../step.hpp"
#include "../shell.hpp"
#include "../print.hpp"
#include "../burst.hpp"
#include "../bruteForce.hpp"
#include "../checks.hpp"
#include "../run_hybGF_P_proj.hpp"
#include "../run_hybGF_PQ_proj.hpp"
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

	int nProj=5;
	double Tsim=1000;
	int Nsamples=500;

//	std::stringstream convert_nProj (argv[1]);
//  	std::stringstream convert_Tsim (argv[2]);
//  	std::stringstream convert_Nsamples (argv[3]

//	if (!(convert_nProj >> nProj ))
//	    exit (EXIT_FAILURE);
//	if (!(convert_Tsim >> Tsim ))
//	    exit (EXIT_FAILURE);
//	if (!(convert_Nsamples >> Nsamples ))
//	    exit (EXIT_FAILURE);

	int stat[3];
	double diffStat[nProj][10];
	double Diff_P [N][Nsamples][nProj];
	double Diff_PQ [N][Nsamples][nProj];
    double Diff_BM [N][Nsamples][nProj];

	//"double diffStat[][10]" where 10 indicates the number of particles
	if (N != 10){
		std::cout << "Error in run_hybGF_P: the number of particles must be equal to 10";
		exit (EXIT_FAILURE);
	}


	for ( int count = 0; count < Nsamples; count++){

        std::cout << count << std::endl;

		for (int d=0; d<3; d++ )
			stat[d] = 0;
	
		for ( int n=0; n<N; n++ )
			for ( int t=0; t<nProj; t++)
				diffStat[t][n] = 0;

		run_hybGF_P_proj ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, nProj, L, stat, diffStat );

		for ( int t=0; t<nProj; t++){

			for ( int n=0; n<N; n++){

				Diff_P [n][count][t] = diffStat[t][n];
	         }

		}


		for (int d=0; d<3; d++ )
			stat[d] = 0;
		for ( int n=0; n<N; n++ )
			for ( int t=0; t<nProj; t++)
				diffStat[t][n] = 0;
	
		run_hybGF_PQ_proj ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, nProj, L, stat, diffStat );

		for ( int n=0; n<N; n++){
			for ( int t=0; t<nProj; t++){

				Diff_PQ [n][count][t] = diffStat[t][n];
		
			}
		}



        for (int d=0; d<3; d++ )
            stat[d] = 0;
        for ( int n=0; n<N; n++ )
            for ( int t=0; t<nProj; t++)
                diffStat[t][n] = 0;

        run_BM ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, Tsim, nProj, L, diffStat );

        for ( int n=0; n<N; n++){
            for ( int t=0; t<nProj; t++){

                Diff_BM [n][count][t] = diffStat[t][n];

            }
        }



    }

	

	double avDiff_P [nProj];
	double avDiff_PQ [nProj];
    double avDiff_BM [nProj];

	for ( int t=0; t<nProj; t++ ){

		avDiff_P [t] = 0;
		avDiff_PQ [t] = 0;
        avDiff_BM [t] = 0;

		for ( int count=0; count<Nsamples; count++){
			for ( int n=0; n<N; n++) {

				avDiff_P [t] += Diff_P [n][count][t]; 
				avDiff_PQ [t] += Diff_PQ [n][count][t];
                avDiff_BM [t] += Diff_BM [n][count][t];

            }
		}

		avDiff_P [t] = avDiff_P [t] / (Nsamples*N);
		avDiff_PQ [t] = avDiff_PQ [t] / (Nsamples*N);
        avDiff_BM [t] = avDiff_BM [t] / (Nsamples*N);

	}


	double sdDiff_P [nProj];
	double sdDiff_PQ [nProj];
    double sdDiff_BM [nProj];

	for ( int t=0; t<nProj; t++ ){

		sdDiff_P [t] = 0;
		sdDiff_PQ [t] = 0;
        sdDiff_BM [t] = 0;

		for ( int count=0; count<Nsamples; count++){
			for (int n=0; n<N; n++){

				sdDiff_P [t] += pow(Diff_P[n][count][t]-avDiff_P[t],2); 
				sdDiff_PQ [t] += pow(Diff_PQ[n][count][t]-avDiff_PQ[t],2);
                sdDiff_BM [t] += pow(Diff_BM[n][count][t]-avDiff_BM[t],2);

            }
		}
	
		sdDiff_P [t] = sqrt(sdDiff_P [t] / pow(Nsamples*N,2) );
		sdDiff_PQ [t] = sqrt(sdDiff_PQ [t] / pow(Nsamples*N,2) );
        sdDiff_BM [t] = sqrt(sdDiff_BM [t] / pow(Nsamples*N,2) );

	}



	std::cout << std::setprecision (7);

	for ( int t=0; t<nProj; t++){

		std::cout << (t+1)*Tsim/nProj << "\t" << avDiff_P[t] << "\t" << avDiff_PQ[t] << "\t" << avDiff_BM[t] << "\t" ;
		std::cout << sdDiff_P[t] << "\t" << sdDiff_PQ[t] << "\t" << sdDiff_BM[t] << "\t" <<(t+1)*Tsim/nProj*6*D_A  << std::endl;

	}

}
