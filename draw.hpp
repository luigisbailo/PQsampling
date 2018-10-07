#pragma once

#include "greensFunct.hpp"

double drawTimeNewt ( double b, double D, double xi ) {

  double t,tmem;

  t = 0.0917517*b*b/D; 
  tmem = -1;

  double S = 1-Sfunct (t,b,D);
  double dS = Sder (t,b,D);
  int count = 0;


  while ( abs(t-tmem) > DRAW_CONVERGENCE | abs(S-xi) > DRAW_CONVERGENCE ) {

    count++;
    if (count > MAX_ITERATIONS){
      std:: setprecision (15);
      std::cout << "Error: finding the root of S was not possible for:" << std::endl;
      std::cout << "b = " << b << "\t D = " << D << "\t xi = " << xi << std::endl;
      return t;
      // exit (EXIT_FAILURE);
    }

    tmem = t;
    t = t + (S-xi)/dS;

    S = 1-Sfunct (t,b,D);
    if ( S==1 ) return t;
    dS = Sder (t,b,D);

  }


  // return count;
  return t;
    
}



double drawPosNewt ( double t, double b, double D, double xi ) {


  double r;

  double t0=0.063;
  double t1=0.234;
  if (t<t0*b*b/D) r = sqrt(t*D)*2;
  else if (t<t1*b*b/D) {
    double R0=2*sqrt(t0*b*b);
    double R1=0.646*b;
    double beta = (R0*exp(pow(t0,0.5)+pow(t1,0.5))-R1*exp(pow(t0,0.5)+pow(t1,0.5)))/(exp(pow(t0,0.5))-exp(pow(t1,0.5)));
    double gamma = -(R0*(exp(pow(t1,0.5))-1)*exp(pow(t0,0.5))-R1*(exp(pow(t0,0.5))-1)*exp(pow(t1,0.5)))  / (exp(pow(t0,0.5))-exp(pow(t1,0.5)));
    r=beta*(1-exp(-pow(t*D/b/b,0.5)))+gamma;
  }
  else r = 0.646*b;


  if (D*t/b/b<0.01){
    std::cout << "ERROR: time too small in drawPos" << std::endl;
  }


  // cout << D*t/b/b << endl;

    double S = Sfunct (t,b,D);
    double P = Pfunct (r,t,b,D,S);
    double dP = Pder (r,t,b,D,S);
    int count = 0;
    double rMem=-1;

 
    while ( abs(r-rMem) > DRAW_CONVERGENCE | abs(P-xi) > DRAW_CONVERGENCE ) {
 

      count++;
      if (count > MAX_ITERATIONS){
        // cout << setprecision (15);
        std::cout << "Error: finding the root of P was not possible for:" << std::endl;
        std::cout << "t = " << t << "\tb = " << b << "\t D = " << D << "\t xi = " << xi << std::endl;
        // exit (EXIT_FAILURE);
        return r;
      }


      rMem = r;
      r = r - (P-xi)/dP;

      // S = Sfunct (t,b,D);
      P = Pfunct (r,t,b,D,S);
      dP = Pder (r,t,b,D,S);

// std::cout << r <<std::endl;
// std::cout << r << "\t" << P << "\t" << dP  <<std::endl;
    }


  // return count;    
  return r;


}



double drawFree ( double t, double D, double xi ){

  double dr = sqrt(6*D*t)/10000;
  double r = dr/2;
  double P=0;


// std::cout << xi << "\t" << t << "\t" << D <<  std::endl;
  
  while (P<xi){

    P += exp(-r*r/4/D/t)/sqrt(4*M_PI*pow(D*t,3))*r*r*dr;
    r += dr;

// std::cout << r << "\t" << P <<"\t" << exp(-r*r/4/D/t)/sqrt(4*M_PI*pow(D*t,3))*r*r*dr << std::endl;

  }

  return r;

}




double drawPosPQ00bis ( double t, double tau, double b, double D, double xi ) {

    if (t>=tau){
        std::cout << t << "\t" << tau << std::endl;
        std::cout << "error drawPosPQ00bis" << std::endl;
        exit (EXIT_FAILURE);
    }
    double q = Sder (tau,b,D);
    double r0;

    if ( t/tau>0.99) {
        double deltaR = 10 * b * (1-t / tau);
        r0 = b - deltaR;
        double PQ = PQ00funct (r0,t,tau,b,D,q);
//        std :: cout << PQ << " " << xi << std::endl;
        while ( PQ > xi ){
            r0 -= deltaR;
            PQ = PQ00funct (r0,t,tau,b,D,q);;
//            std :: cout << PQ << std::endl;

        }

    }
    else {
        r0 = 0;
      }

    double r1 = b;
    double r = (r0 + r1) /2;

    double P = PQ00funct (r,t,tau,b,D,q);
    int count = 0;
    double rMem = -1;

    while ( abs(P-xi)>DRAW_CONVERGENCE | abs(r-rMem) > DRAW_CONVERGENCE*b ){
//        std::cout << r0 << " " << r1 << " " << r << " " << P << " " << xi << std::endl;

        count++;
        if (count > MAX_ITERATIONS){
            std::cout << std::setprecision (15);
            std::cout << "Error: finding the root of P was not possible for:" << std::endl;
            std::cout << "t = " << t << "\t b = " << b << "\t D = " << D << "\t xi = " << xi << std::endl;
            exit (EXIT_FAILURE);
            }

        rMem = r;
        if ( P<xi ) r0 = r;
        else r1 = r;

        r = (r0 + r1)/2;
        P = PQ00funct (r,t,tau,b,D,q);

    };

    return r;

}


double drawPosPQbis ( double t, double radius0, double tau, double b, double D, double xi ) {

    if (t>=tau){
        std::cout << "error time drawPosPQbis" << std::endl;
        exit (EXIT_FAILURE);
    }
    if (radius0>=b){
        std::cout << "error position drawPosPQbis" << std::endl;
        exit (EXIT_FAILURE);
    }
    double r0,r1,r,q,P;
    r = radius0;

    q = qFunct (radius0,tau,b,D);
    P = PQfunct (r,t,radius0,tau,b,D,q);
//    std::cout<<"\t" << xi << "\t" << P << "\t" << P-xi <<std::endl;

//    double deltaR = (b - radius0)/4;

    if (P>xi){
        r0 = 0;
        r1 = r;
    }
    else{
        r0 = r;
        r1 = b;
    }

    r = (r0 + r1) /2;
    P = PQfunct (r,t,radius0,tau,b,D,q);
//    std::cout << r0 << "\t" <<r  << "\t" << r1 << "\t" << P-xi << std::endl;

    int count = 0;
    double rMem = -1;

    while ( abs(r-rMem) >DRAW_CONVERGENCE*b && abs(P-xi) > DRAW_CONVERGENCE ){

        count++;
        if (count > MAX_ITERATIONS){
            std::cout << std::setprecision (15);
            std::cout << "Error: finding the root of PQ was not possible for:" << std::endl;
            std::cout << "r0 = " << radius0 << "\tt = " << t << "\ttau = " << tau << "\t b = " << b << "\t D = " << D << "\t xi = " << xi << std::endl;
//             exit (EXIT_FAILURE);
            return r;
        }

        rMem = r;

        if ( P<xi ) r0 = r;
        else r1 = r;

        r = (r0 + r1)/2;
        P = PQfunct (r,t,radius0,tau,b,D,q);

//         std::cout << r0 << "\t" <<r  << "\t" << r1 << "\t" << P-xi << std::endl;

    };

    return r;

}


