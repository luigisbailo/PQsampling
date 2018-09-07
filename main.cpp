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
#include <random>

using namespace std::chrono;

#include "parameters.hpp"
#include "tools.hpp"
#include "greensFunct.hpp"
#include "draw.hpp"
#include "init.hpp"
#include "step.hpp"
#include "shell.hpp"
#include "print.hpp"
#include "burst.hpp"
#include "bruteForce.hpp"
#include "checks.hpp"

int main () {         

// cout << drawTimeNewt (1.5541445463895798e-08,0.001,0.38271502428688109) << endl;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  high_resolution_clock::time_point t2,t1;

  double tau = 20;
  double t = 0.5*tau;
  double b = 1.54;
  double D = 0.01;
  double radius0 = 0.1;

  double dX = 0.01;

//    for ( double X = dX; X<1; X += dX){
//
////      std::cout << X << "\t" << drawPosPQ00Newt (t, tau, b, D, X) << "\t" << drawPosPQ00Newt (t, tau, b, D, X) << std::endl;
//        std::cout << X << "\t" << PQ00funct (X, t, tau, b, D, Sder(t,b,D)) << std::endl;
//    }

//    int nsamples = 10000;
//
//  t1 = high_resolution_clock::now();
//  for (int count = 0; count<nsamples; count++ ){
//      double xi =  distribution(generator);
//      double r = drawPosNewt(t,b,D,xi);
////      std::cout << r << std::endl;
//    }
//  t2 = high_resolution_clock::now();
//
//  double t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;
//
//  std::cout << "P\t" << t12 << std::endl;
//
//  t1 = high_resolution_clock::now();
//  for (int count = 0; count<nsamples; count++ ){
//    double xi =  distribution(generator);
//    double r = drawPosPQ00bis(t,tau,b,D,xi);
////      std::cout << r << std::endl;
//  }
//  t2 = high_resolution_clock::now();
//
//  t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;
//
//  std:: cout << "PQ bis\t" << t12 << std::endl;
//


  double dr = 0.01;
  double R;
  double PQ_temp;
  b = 2;
  radius0 = 1.9;
  double q = qFunct (radius0,tau,b,D);
  t=0.7*tau;
  int nsamples = 10000;
//  drawPosPQbis(t,radius0, tau,b,D,0.1896897718996439219374395);
  t1 = high_resolution_clock::now();
  for (int count = 0; count<nsamples; count++ ){
      double xi =  distribution(generator);
      double r = drawPosPQbis(t,radius0, tau,b,D,xi);
//      std::cout << std::endl << std::endl;
    }
  t2 = high_resolution_clock::now();

  double t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

  std::cout << "PQ\t" << t12 << std::endl;

  t1 = high_resolution_clock::now();
  for (int count = 0; count<nsamples; count++ ){
      double xi =  distribution(generator);
      double r = drawPosNewt(t,b,D,xi);
//      std::cout << r << std::endl;
    }
  t2 = high_resolution_clock::now();

  t12 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count())/1000000;

  std::cout << "P\t" << t12 << std::endl;

//  for (double r=dr; r<b; r+=dr) {
//    std::cout << r << "\t" << PQfunct(r, t, radius0, tau, b, D, q)<<std::endl;
//  }

//  double PQ_max=0;
//  double R_max=0;
//


//  for (dX = 0; dX<1; dX += 0.001) {
//    std::cout<< dX << "\t" << drawPosPQbis(t, radius0, tau, b, D, dX) << std::endl << std::endl;
//  }
//
//  for ( double r=dr; r<b; r+=dr ) {
//
////      PQ_temp = PQ00funct(r,t,tau,b,D,Sder(tau,b,D));
//
//    PQ_temp = PQ00der(r,t,tau,b,D,Sder(tau,b,D));
//    if ( PQ_temp > PQ_max ){
//      R_max = r;
//      PQ_max = PQ_temp;
//    }
////    std :: cout << r << "\t" << PQ_temp << std::endl;
//  }
//
//  std::cout << R_max << std::endl;
//  std::cout << t/tau*b << std::endl;

}




