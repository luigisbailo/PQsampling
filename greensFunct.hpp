// author luigisbailo


#pragma once

double Sfunct ( double t, double b, double D) {
//This function returns 1-Survival probability, i.e. the probability of having already left the domain

  double coeff1 = exp ( -M_PI*M_PI*D*t/(b*b) );
  double coeff2 = 2;
  double S = 0;
  double term,termA,termB;
  int m = 1;
  double conv = 0.0000001/(b*b/D);

  if ( t>= b*b/D/100 ) {

    do {

      termA = pow (coeff1,m*m);
      termB = pow (coeff1,(m+1)*(m+1));
      term = coeff2 * (termA-termB);
      S += term;
      m += 2;

    } while ( abs (term) > conv );

    S = 1-S;

  }
  else 
    S = 0;

  return S;

}


double Sder ( double t, double b, double D) {

  double coeff1 = exp ( - M_PI*M_PI*D*t/(b*b) );
  double coeff2 = 2*D*M_PI*M_PI/b/b;
  double S;
  double term;
  double termA,termB,termC,termD;
  int m;

  double conv = 0.00001/(b*b/D);

  S = 0;
  m = 1;

  if ( t>= b*b/D/100 ) {

    do {

      termA = pow (coeff1,m*m);
      termB = m*m;
      termC = pow (coeff1,(m+1)*(m+1));
      termD = (m+1)*(m+1);
      term = coeff2 * ( termA*termB - termC*termD );

      S += term;
      m += 2;

    } while ( abs (term) > conv | termC > termD/100 | S < 0 );
    //The second condition takes into account that the series, for short times, can be negative in the beginning, because of the term m*m 
    //The series, when the term m*m  is dominant at short times, grows to a peak, and then finally deceases to zero
    //The S<0 condition is to ensure a longer convergence for short times
  }
  else 
    S = 0;

  return S;

}


double Pfunct ( double radius, double t, double b, double D, double S ) {

  double coeff1 = exp ( - M_PI*M_PI*D*t/(b*b));
  double coeff2 = 2/(b*M_PI);
  double P;
  double  term, termA, termB;
  int m;
  double conv = 0.00000001/b;

  P = 0;
  m = 1;

  if (t < b*b/D*0.2) {      

    do {

      termA =  pow(coeff1,m*m);
      termB =  b*sin(m*M_PI*radius/b)/m - M_PI*radius*cos(m*M_PI*radius/b); 
      term =  coeff2 * termA * termB / (1-S);
      P += term;
      m += 1;

    } while ( abs(term) > conv  |  termA>termB/100 |  P>1 );

  } else {

    P = sin(M_PI*radius/b)/M_PI - cos(M_PI*radius/b)*radius/b;

  }

  return P;

}


double Pder ( double radius, double t, double b, double D, double S ) {

  double coeff1 = exp ( -M_PI*M_PI*D*t/(b*b)  );
  double coeff2 = 2*M_PI/b/b;
  double P;
  double term;
  double termA,termB;
  int m;
  double conv = 0.0001/b;

  P = 0;
  m = 1;

  if (t < b*b/D*0.2) {      

    do {

      termA = pow(coeff1,m*m);
      termB = sin (m*M_PI*radius/b) * m * radius;
      term = coeff2 * termA * termB / (1-S);
      P += term;
      m += 1;

    } while ( abs(term) > conv  |  termA>termB/100 );

  }
  else {

    P = M_PI*sin(M_PI*radius/b)*radius/b/b;

  }

  return P;

}


double PQfunct ( double radius, double t, double r0, double tau, double b, double D, double q ) {

  double coeff1 = exp ( - M_PI*M_PI*D/(b*b));
  double coeff2 = 2*D/b/r0;
  double PQ;
  double  term1, term2, termA, termB;
  int m,n;
  double conv = 0.00000001/b;


  PQ = 0;
  m = 1;
      
  do {

    n = 1;
    term2 = 0;

    do {

      termA =  pow( coeff1, m*m*t + n*n*(tau-t) );

      if ( m != n)
        termB =  n*sin((m-n)*M_PI*radius/b)/(m-n) - n*sin((m+n)*M_PI*radius/b)/(m+n);       
      else if ( m == n )
        termB =  n*M_PI*radius/b - n*sin((m+n)*M_PI*radius/b)/(m+n); 
      termB = termB * sin(m*M_PI*r0/b);
      term1 =  coeff2 * termA * termB / q;
      term2 += term1;
      n += 1;
      termA =  pow( coeff1, m*m*t + n*n*(tau-t) );

      if ( m != n)
        termB =  n*sin((m-n)*M_PI*radius/b)/(m-n) - n*sin((m+n)*M_PI*radius/b)/(m+n);       
      else if ( m == n )
        termB =  n*M_PI*radius/b - n*sin((m+n)*M_PI*radius/b)/(m+n); 
      termB = termB * sin(m*M_PI*r0/b);

      term1 =  coeff2 * termA * termB / q;
      term2 -= term1;
      n += 1;

    } while (  abs (term1) > conv | termA*n*n>0.000001 | n<10 ) ;

    PQ += term2;
    m += 1;

  } while ( abs (term2) > conv | termA*m*m>0.00001 | m<10 );

  return PQ;

}


double PQder ( double radius, double t, double r0, double tau, double b, double D, double q ) {

  double coeff1 = exp ( -M_PI*M_PI*D/(b*b) );
  double coeff2 = 2*M_PI*D/b/b/r0;
  double PQ;
  double  term1, term2, termA, termB;
  int m,n;
  double conv = 0.00000001/b;

  PQ = 0;
  m = 1;

  do {

    n = 1;
    term2 = 0;

    do {

      termA = pow (coeff1, m*m*t + n*n*(tau-t) );

      termB = n*cos((m-n)*M_PI*radius/b) - n*cos((m+n)*M_PI*radius/b);
      termB = termB * sin(m*M_PI*r0/b);

      n += 1;

      term1 =  coeff2 * termA * termB / q;
      term2 += term1;


      termA = pow (coeff1, m*m*t + n*n*(tau-t) );

      termB = n*cos((m-n)*M_PI*radius/b) - n*cos((m+n)*M_PI*radius/b);
      termB = termB * sin(m*M_PI*r0/b);

      n += 1;

      term1 =  coeff2 * termA * termB / q;
      term2 -= term1;

    } while ( abs(term1) > conv | termA*n*n>0.0000001 | n<200 );

    PQ += term2;
    m += 1;

  } while ( abs(term2) > conv | termA*m*m>0.00001 | m<2000);

  return PQ;

}


double PQ00funct ( double radius, double t, double tau, double b, double D, double q ) {


  double coeff1 = exp ( - M_PI*M_PI*D/(b*b));
  double coeff2 = 2*M_PI*D/b/b;
  double PQ;
  double  term1, term2, termA, termB;
  int m,n;
  double conv = 0.00000001/b;

 PQ = 0;
  m = 1;
      
  do {

    n = 1;
    term2 = 0;

    do {

      termA =  pow( coeff1, m*m*t + n*n*(tau-t) );

      if ( m != n)
        termB =  m*n*sin((m-n)*M_PI*radius/b)/(m-n) - m*n*sin((m+n)*M_PI*radius/b)/(m+n);       
      else if ( m == n )
        termB =  m*n*M_PI*radius/b - m*n*sin((m+n)*M_PI*radius/b)/(m+n); 

      term1 =  coeff2 * termA * termB / q;
      term2 += term1;
      n += 1;
      termA =  pow( coeff1, m*m*t + n*n*(tau-t) );

      if ( m != n)
        termB =  m*n*sin((m-n)*M_PI*radius/b)/(m-n) - m*n*sin((m+n)*M_PI*radius/b)/(m+n);       
      else if ( m == n )
        termB =  m*n*M_PI*radius/b - m*n*sin((m+n)*M_PI*radius/b)/(m+n); 

      term1 =  coeff2 * termA * termB / q;
      term2 -= term1;

      n += 1;

    } while (  abs (term1) > conv | termA*n*n>0.000001 | n<5 ) ;

    PQ += term2;
    m += 1;

  } while ( abs (term2) > conv | termA*m*m>0.000001 | m<5 );
 
  return PQ;

}


double PQ00der ( double radius, double t, double tau, double b, double D, double q ) {

  double coeff1 = exp ( -M_PI*M_PI*D/(b*b) );
  double coeff2 = 2*M_PI*M_PI*D/b/b/b;
  double PQ;
  double  term1, term2, termA, termB;
  int m,n;
  double conv = 0.00000001/b;

  PQ = 0;
  m = 1;

  do {

    n = 1;
    term2 = 0;

    do {

      termA = pow (coeff1, m*m*t + n*n*(tau-t) );

      termB = m*n*cos((m-n)*M_PI*radius/b) - m*n*cos((m+n)*M_PI*radius/b);

      n += 1;

      term1 =  coeff2 * termA * termB / q;
      term2 += term1;


      termA = pow (coeff1, m*m*t + n*n*(tau-t) );

      termB = m*n*cos((m-n)*M_PI*radius/b) - m*n*cos((m+n)*M_PI*radius/b);

      n += 1;

      term1 =  coeff2 * termA * termB / q;
      term2 -= term1;

    } while ( abs(term1) > conv | termA*n*n>0.000001 | n<5 );

    PQ += term2;
    m += 1;
  } while ( abs(term2) > conv | termA*m*m>0.000001 | m<15);

  return PQ;

}


double qFunct ( double radius, double t, double b, double D ) {

  double coeff1 = exp ( - M_PI*M_PI*D*t/(b*b));
  double coeff2 = 2*M_PI*D/radius/b;
  double q;
  double  term, term1, term2;
  int m;
  double conv = 0.000000001/(b*b/D);

  q = 0;
  m = 1;
      
  do {

    term1 =  pow(coeff1,m*m)*sin(m*M_PI*radius/b)*m;
    term2 = pow(coeff1,(m+1)*(m+1))*sin((m+1)*M_PI*radius/b)*(m+1); 
    term =  coeff2 * (term1 - term2);
    q += term;
    m += 2;


  } while ( abs(term) > conv | m<100 );

  return q;

}


double pFunct ( double radius, double t, double b, double D ) {

  double coeff1 = exp ( -M_PI*M_PI*D*t/(b*b)  );
  double coeff2 = 1./2/b/b;
  double P;
  double term;
  double termA,termB;
  int m;
  double conv = 0.0000001/b;

  P = 0;
  m = 1;

  do {

    termA = pow(coeff1,m*m);
    termB = sin (m*M_PI*radius/b) * m / radius;
    term = coeff2 * termA * termB ;
    P += term;
    m += 1;

  } while ( abs(term) > conv  |  m<20 | P<0);

  return P;

}


double pFunct_init ( double r, double t, double r0, double t0, double b, double D ) {

  double coeff1 = exp ( -M_PI*M_PI*D*(t-t0)/(b*b)  );
  double coeff2 = 1./2/M_PI/b;
  double P;
  double term;
  double termA,termB;
  int m;
  double conv = 0.0000001/b;

  P = 0;
  m = 1;

  do {

    termA = pow(coeff1,m*m);
    termB = sin (m*M_PI*r/b) * sin (m*M_PI*r0/b)  / r / r0;
    term = coeff2 * termA * termB ;
    P += term;
    m += 1;

  } while ( abs(term) > conv  |  m<20 );

  return P;

}


double Sfunct_init ( double tau, double radius, double t, double b, double D) {
//This function returns 1-Survival probability, i.e. the probability of having already left the domain

  double coeff1 = exp ( -M_PI*M_PI*D*(tau-t)/(b*b) );
  double coeff2 = 2*b/M_PI;
  double S = 0;
  double term,termA,termB;
  int m = 1;
  double conv = 0.0000001/(b*b/D);

  if ( tau-t>= b*b/D/100 ) {

    do {

      termA = pow (coeff1,m*m)*sin(m*M_PI*radius/b)/radius/m;
      termB = pow (coeff1,(m+1)*(m+1))*sin((m+1)*M_PI*radius/b)/radius/(m+1);
      term = coeff2 * (termA-termB);


      S += term;
      m += 2;

    } while ( abs (term) > conv );

    S = 1-S;
  }
  else 
    S = 0;

  return S;

}

double Sder2 ( double t, double b, double D) {

  double coeff1 = exp ( - M_PI*M_PI*D*t/(b*b) );
  double coeff2 = -2*D*D*M_PI*M_PI*M_PI*M_PI/(b*b*b*b);
  double S;
  double term;
  double termA,termB,termC,termD;
  int m;

  double conv = 0.0001/(b*b/D);

  S = 0;
  m = 1;

    do {

        termA = pow (coeff1,m*m);
        termB = m*m*m*m;
        termC = pow (coeff1,(m+1)*(m+1));
        termD = (m+1)*(m+1)*(m+1)*(m+1);
        term =  ( termA*termB - termC*termD );

        S += term;
        m += 2;

    } while ( abs (term) > conv | termC > termD/100  );

  return S*coeff2;

}


double Pder2 ( double radius, double t, double b, double D, double S ) {

  double coeff1 = exp ( -M_PI*M_PI*D*t/(b*b )  );
  double coeff2 = 2*M_PI/b/b/(1-S);
  double P;
  double termA, termB, term;
  int m;
  double conv = 0.0001/b;

  term = 0;
  P = 0;
  m = 1;

  do {

    termA = pow(coeff1,m*m);
    termB =  m * sin (m*M_PI*radius/b) + m*m*M_PI*radius/b * cos (m*M_PI*radius/b);
    term  = termA*termB;
    P += term;
    m += 1;

  } while ( abs (term) > conv | termA>termB/100 );

  return P*coeff2;

}

