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
        exit (EXIT_FAILURE);
        std::cout << "error drawPosPQ00bis" << std::endl;
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

    double r0 = 0;
    double r1 = b;
    double r = radius0;
    double q = qFunct (radius0,tau,b,D);
    double P = PQfunct (r,t,radius0,tau,b,D,q);
    // double dP = PQder (r,t,radius0,tau,b,D,q);
    // double rMem=-1;
    int count = 0;

    // if ( xi<pow(10,-16) )
    //   return 0;

    while ( abs(P-xi)>DRAW_CONVERGENCE*b | abs(P-xi) > DRAW_CONVERGENCE ){

        count++;
        if (count > MAX_ITERATIONS){
            std::cout << std::setprecision (15);
            std::cout << "Error: finding the root of P was not possible for:" << std::endl;
            std::cout << "r0 = " << radius0 << "\tt = " << t << "\ttau = " << tau << "\t b = " << b << "\t D = " << D << "\t xi = " << xi << std::endl;
            // exit (EXIT_FAILURE);
            return r;
        }

        if ( P<xi ) r0 = r;
        else r1 = r;

        r = (r0 + r1)/2;
        P = PQfunct (r,t,radius0,tau,b,D,q);



        // rMem =r;
        // r = r - (P-xi)/dP;

        // q = qFunct (radius0,tau,b,D);
        // P = PQfunct (r,t,radius0,tau,b,D,q);
        // dP = PQder (r,t,radius0,tau,b,D,q);



        // std::cout << r0 << "\t" <<r  << "\t" << r1 << "\t" << P-xi << std::endl;


    };

    return r;

}


//double drawPosPQ00Newt ( double t, double tau, double b, double D, double xi ) {
//
//
//    double r= b*t/tau;
//
//    if (D*t/b/b<0.01){
//        std::cout << "ERROR: time too small in drawPos" << std::endl;
//    }
//
//    double q = Sder (tau,b,D);
//    double P = PQ00funct (r,t,tau,b,D,q);
//    double dP = PQ00der (r,t,tau,b,D,q);
//    int count = 0;
//    double rMem=-1;
//
//
//    while ( abs(r-rMem) > DRAW_CONVERGENCE*b | abs(P-xi) > DRAW_CONVERGENCE ) {
//
//
//        count++;
//        if (count > MAX_ITERATIONS){
//            // cout << setprecision (15);
//            std::cout << "Error: finding the root of P was not possible for:" << std::endl;
//            std::cout << "t = " << t << "\tb = " << b << "\t D = " << D << "\t xi = " << xi << std::endl;
//            // exit (EXIT_FAILURE);
//            return r;
//        }
//
//
//        rMem = r;
//        r = r - (P-xi)/dP;
//
//        P = PQ00funct (r,t,tau,b,D,q);
//        dP = PQ00der (r,t,tau,b,D,q);
//
//// std::cout << r <<std::endl;
//// std::cout << r << "\t" << P << "\t" << dP  <<std::endl;
//    }
//
//
//    // return count;
//    return r;
//
//
//}
//




// double drawFree ( double t, double D, double xi){

//   double dr = sqrt(t*D)/100;
//   double R=0;
//   double r = dr/2;

//   while (R<xi){

//     R += r*r*exp(-r*r/4/D/t)/sqrt(4*M_PI*pow(D*t,3))*dr;
//     r += dr;

//   }

//   return R;

// }



// double draw0PosDer2Bis ( double t, double b, double D ) {

//   double S = Sfunct(t,b,D);
//   double xi = 0.001;

//     while ( Pder2 (xi*b,t,b,D,S) >0){ 
//     xi += 0.001;
//     }

//   // std::cout << xi*b << "\t" << Pder2 (xi*b,t,b,D,S) << std::endl;
//   // if (t<0.25*b*b/D) std::cout << 2.584*D*t/b << std::endl;
//   // else std::cout << 0.646*b << std::endl;

//   return xi*b;

// }


// double drawPosBis ( double t, double b, double D, double xi ) {

//   double r,P;
//   double S = Sfunct(t,b,D);
//   double r0 = 0;
//   double r1 = b;


//   int count =0;
//   do{

//   count++;
//     r = (r0 + r1)/2; 
//     P = Pfunct (r,t,b,D,S);
//   std::cout << r0 << "\t" << r << "\t" << r1 << "\t" << P-xi << std::endl;

//     if ( P<xi ) r0 = r;
//     else r1 = r;  

//   } while ( abs(P-xi)>DRAW_CONVERGENCE && count<30); 

//   // return count;

//   return r;

    
// }





// double drawTime ( double b, double D, double xi, SfunctVar *mySvar) {

//   struct Sinters_params params = { b, D, xi, mySvar};

//   (*mySvar).interval = (b*b/D)/20;
//   //the interval size is given considering the average time exit b^2/(D*6*pi)

//   (*mySvar).F.function = &Sinters;
//   (*mySvar).F.params = &params;
    
//   (*mySvar).solver_type = gsl_root_fsolver_brent;
//   (*mySvar).solver = gsl_root_fsolver_alloc((*mySvar).solver_type);
    
//   (*mySvar).Sinters_low=0; 
//   (*mySvar).Sinters_high=0.+(*mySvar).interval;  

//   (*mySvar).status = GSL_CONTINUE;

  
//   for (int i = 0; i <= MAX_ITERATIONS && (*mySvar).status == GSL_CONTINUE; i++){

//     if ( (GSL_FN_EVAL(&(*mySvar).F,(*mySvar).Sinters_low) * GSL_FN_EVAL(&(*mySvar).F,(*mySvar).Sinters_high)) < 0 ) {
	
//       gsl_root_fsolver_set ((*mySvar).solver, &(*mySvar).F, (*mySvar).Sinters_low, (*mySvar).Sinters_high);
//       (*mySvar).status = GSL_CONTINUE;
//       for (int i=1; 1<100 && (*mySvar).status==GSL_CONTINUE; i++) {
  
//       	(*mySvar).status = gsl_root_fsolver_iterate ((*mySvar).solver);
//       	if ((*mySvar).status != GSL_SUCCESS)  break;
//       	//!! TO EVALUATE THE POSITION OF THIS BREAK	  
//       	(*mySvar).r = gsl_root_fsolver_root ((*mySvar).solver);
//       	(*mySvar).Sinters_low = gsl_root_fsolver_x_lower ((*mySvar).solver);
//       	(*mySvar).Sinters_high = gsl_root_fsolver_x_upper ((*mySvar).solver);
//       	(*mySvar).status = gsl_root_test_interval ((*mySvar).Sinters_low, (*mySvar).Sinters_high, 0, 0.001);
	  
//       }
//     }

//     (*mySvar).Sinters_high = (*mySvar).Sinters_low + (*mySvar).interval;      
//     (*mySvar).Sinters_low+=(*mySvar).interval; 
//     (*mySvar).Sinters_high+=(*mySvar).interval;
      
//   }

//   if ((*mySvar).status == GSL_CONTINUE) {
//     printf("ERROR: too many iterations in drawTime\n");
//     exit (EXIT_FAILURE);
//   }
  
//     gsl_root_fsolver_free ((*mySvar).solver);
//     return  (*mySvar).r;

// };


// double drawPos_P ( double t, double b, double D,  double xi, SfunctVar *mySvar, PfunctVar *myPvar ) {

//   if (t*D<0.0000000001) return 0;
//   if (1-xi<0.000000000001) return b;


//   struct P_inters_params params = { t, b, D, xi, mySvar, myPvar};

//   (*myPvar).interval = b/20;
//   //TO VERIFY THE INTERVAL

//   (*myPvar).F.function = &P_inters;
//   (*myPvar).F.params = &params;
    
//   (*myPvar).solver_type = gsl_root_fsolver_brent;
//   (*myPvar).solver = gsl_root_fsolver_alloc((*myPvar).solver_type);
    
//   (*myPvar).Pinters_low=0; 
//   (*myPvar).Pinters_high=0.+(*myPvar).interval;  

//   (*myPvar).status = GSL_CONTINUE;


//   for (int i = 0; i <= MAX_ITERATIONS && (*myPvar).status == GSL_CONTINUE; i++){
      
//       if ( (GSL_FN_EVAL(&(*myPvar).F,(*myPvar).Pinters_low) * GSL_FN_EVAL(&(*myPvar).F,(*myPvar).Pinters_high)) < 0 ) {
	
//       	gsl_root_fsolver_set ((*myPvar).solver, &(*myPvar).F, (*myPvar).Pinters_low, (*myPvar).Pinters_high);
//       	(*myPvar).status = GSL_CONTINUE;
//       	for (int i=1; 1<100 && (*myPvar).status==GSL_CONTINUE; i++) {
//     	  (*myPvar).status = gsl_root_fsolver_iterate ((*myPvar).solver);
//     	  if ((*myPvar).status != GSL_SUCCESS)  break;
//     	  //!! TO EVALUATE THE POSITION OF THIS BREAK	  
//     	  (*myPvar).r = gsl_root_fsolver_root ((*myPvar).solver);
//     	  (*myPvar).Pinters_low = gsl_root_fsolver_x_lower ((*myPvar).solver);
//     	  (*myPvar).Pinters_high = gsl_root_fsolver_x_upper ((*myPvar).solver);
//     	  (*myPvar).status = gsl_root_test_interval ((*myPvar).Pinters_low, (*myPvar).Pinters_high, 0, 0.001);
    	  
//   	  }
//     }

//   (*myPvar).Pinters_high = (*myPvar).Pinters_low + (*myPvar).interval;      
//   (*myPvar).Pinters_low += (*myPvar).interval; 
//   (*myPvar).Pinters_high += (*myPvar).interval;
  
//   } while ((*myPvar).status == GSL_CONTINUE);

//   if ((*myPvar).status == GSL_CONTINUE) {
//     printf("ERROR: too many iterations in drawPos_P\n");
//     exit (EXIT_FAILURE);
//   }


    
//     gsl_root_fsolver_free ((*myPvar).solver);
//     return (*myPvar).r;

// };


