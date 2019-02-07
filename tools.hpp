#pragma once

#include "draw.hpp"

class particle {
public:
    int label;
    int type;
    double pos [3];
    double pos_init [3];
    int pos_period [3];
    double shell;
    double time;
    double tau_exit;
    double tau_exitSampled;
    double pos_exit [3];
    double radius;
    double Diff;
    double sqrtDiff;
    double R_bd;
    double R_gfrd;
    double burstR;
    double shellRed;
    bool burst;
    bool gf;
    double distBM;
    double distGF;
    double distBurst;
    double displPQ [3][MEMORY_PQ];
    int countPQ;
    int totPQdispl;
    bool active;

 };


  class BFdistances {
  public: 
    double dd[10];

  };


double dist2_per (particle *A, particle *B, double L ) {
//squared particles distance in periodic boundary conditions
 
 double XYZ [3];

  XYZ[0] = A->pos[0]-B->pos[0];
  if ( abs(XYZ[0]) > L/2 ) XYZ[0] = L-abs(XYZ[0]);
  XYZ[1] = A->pos[1]-B->pos[1];
  if ( abs(XYZ[1]) > L/2 ) XYZ[1] = L-abs(XYZ[1]);
  XYZ[2] = A->pos[2]-B->pos[2];
  if ( abs(XYZ[2]) > L/2 ) XYZ[2] = L-abs(XYZ[2]);

  return XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2];  

}

double dist2next_per (particle *A, particle *B, double L ) {

  double XYZ [3];

  XYZ[0] = A->pos[0]-B->pos_exit[0];
  if ( abs(XYZ[0]) > L/2 ) XYZ[0] = L-abs(XYZ[0]);
  XYZ[1] = A->pos[1]-B->pos_exit[1];
  if ( abs(XYZ[1]) > L/2 ) XYZ[1] = L-abs(XYZ[1]);
  XYZ[2] = A->pos[2]-B->pos_exit[2];
  if ( abs(XYZ[2]) > L/2 ) XYZ[2] = L-abs(XYZ[2]);

  return XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2];

}

void printPos_per (particle *particles, int *partList, int N){

  std::cout <<"Part.\t"<<"Time\t\t"<<"Tau\t\t"<<"Shell\t"
       <<"x\t"<<"y\t"<<"z\t"
       <<"x_ex\t"<<"y_ex\t"<<"z_ex\t"
       <<"x_per\t"<<"y_per\t"<<"z_per\t"
       <<"Radius\t" <<"BURST\t"<<"GF\n";

  for (int count=0; count<N; count++){
    std::cout <<particles[partList[count]].label << "\t"
   <<particles[partList[count]].time << "\t\t" 
   <<particles[partList[count]].tau_exit<<"\t\t"
   <<particles[partList[count]].shell<<"\t"
   <<particles[partList[count]].pos[0]<<"\t"
   <<particles[partList[count]].pos[1]<<"\t"
   <<particles[partList[count]].pos[2] <<"\t"
   <<particles[partList[count]].pos_exit[0]<<"\t"
   <<particles[partList[count]].pos_exit[1]<<"\t"
   <<particles[partList[count]].pos_exit[2] <<"\t"
   <<particles[partList[count]].pos_period[0]<<"\t"
   <<particles[partList[count]].pos_period[1]<<"\t"
   <<particles[partList[count]].pos_period[2] <<"\t"
   <<particles[partList[count]].pos_init[0]<<"\t"
   <<particles[partList[count]].pos_init[1]<<"\t"
   <<particles[partList[count]].pos_init[2] <<"\t"
   <<particles[partList[count]].radius <<"\t"
   <<particles[partList[count]].burst <<"\t"
   <<particles[partList[count]].gf<<"\n";

  }
}


void printPos_annih (particle *particles, int *partList, int N){

  std::cout <<"Part.\t"<<"Time\t\t"<<"Tau\t\t"
            <<"Radius\t" <<"BURST\t"<<"GF\t"<<"Active\n";

  for (int count=0; count<N; count++){
    std::cout <<particles[partList[count]].label << "\t"
              <<particles[partList[count]].time << "\t\t"
              <<particles[partList[count]].tau_exit<<"\t\t"
              <<particles[partList[count]].radius <<"\t"
              <<particles[partList[count]].burst <<"\t"
              <<particles[partList[count]].gf<<"\t"
              <<particles[partList[count]].active << std::endl;
  }
}


void getDist ( particle *particles, int* partList, double *distRow, double *maxSh ,int N, double L ) {

  int pos;
  //The j-th element of the vector "distRow[]" is filled with the distance with the j-th particle
  //The j-th particle means the j-th ordered on a exit-time basis
  //The 0 position is left empty, since "getDist()" is always called for the particle 0   
  for ( int j=1; j<N; j++ ) {

    if (particles[partList[j]].active == false) {
      distRow [j] = L;
    }
    else if (particles[partList[j]].gf){
      distRow [j] = sqrt(dist2_per ( &particles[partList[0]], &particles[partList[j]], L )) - particles[partList[0]].radius - particles[partList[j]].radius;
    }
    //distance in case of BM. The two particles times are in within an integration step
    else if ( abs(particles[partList[0]].time-particles[partList[j]].time) < abs(particles[partList[0]].time-particles[partList[j]].tau_exit) ){
      distRow [j] = sqrt(dist2_per ( &particles[partList[0]], &particles[partList[j]], L )) - particles[partList[0]].radius - particles[partList[j]].radius;
    }  
    else {
      distRow [j] = sqrt(dist2next_per ( &particles[partList[0]], &particles[partList[j]], L )) - particles[partList[0]].radius - particles[partList[j]].radius;
    }

  }

  distRow [0] = distRow[1];
  pos = 1;

  //"distRow[0]" is filled with the lowest distance; "pos" labels the position of the closest particle  
  for (int j=2; j<N; j++){

    if (distRow[j]<distRow[0] ){
      distRow[0] = distRow[j];
      pos = j;    

    }
  }

  // "maxSh" says the largest shell that the particle 0 can construct 
  *maxSh = distRow[pos] / ( 1 + particles[partList[pos]].sqrtDiff/particles[partList[0]].sqrtDiff );


}



bool compareTime (particle A, particle B)  {return A.tau_exit < B.tau_exit;}


void sortPart ( particle *particles, int *partList, int N) {

  int tempList;

  for (int n=0; n<N-1; n++){

    if (particles[partList[n]].tau_exit>particles[partList[n+1]].tau_exit){

      tempList=partList[n];
      partList[n]=partList[n+1];
      partList[n+1]=tempList;

    }
    else break;

  }

}

void printPos_per (particle *particle, int *partList, int N);


void sortBurst ( particle *particles, int *partList, int N) {

  int tempList;

  for (int n=1; n<N; n++){

    if (particles[partList[n]].burst){

      tempList = partList[n];
      for ( int m = n; m>1; m-- ){

        partList[m] = partList[m-1];

      }

      partList[1]=tempList;
    }

  }

}



void checkBound (double *pos, int *pos_period, double L) {

  if ( pos[0] > L ) {
    pos[0] -= L;
    pos_period[0] ++;
  }
  else if ( pos[0] < 0 ) {
    pos[0] += L;
    pos_period[0] --;
  }

  if ( pos[1] > L ) {
    pos[1] -= L;
    pos_period[1] ++;
  }
  else if ( pos[1] < 0 ) {
    pos[1] += L;
    pos_period[1] --;
  }

  if ( pos[2] > L ) {
    pos[2] -= L;
    pos_period[2] ++;
  }
  else if ( pos[2] < 0 ) {
    pos[2] += L;
    pos_period[2] --;
  }

};



void polarTransf ( double *pos, double R, double u, double v ) {
//Theta is defined in [0,2pi]
//Phi is defined in [0,pi]
  
  double theta = 2 * M_PI * u;
  double phi = acos( 2*v - 1 );

  pos[0] = R * cos(theta) * sin(phi);      
  pos[1] = R * sin(theta) * sin(phi);
  pos[2] = R * cos(phi);  

}



void polarTransf_angles ( double *pos, double R, double theta, double phi ) {
//Theta is defined in [0,2pi]
//Phi is defined in [0,pi]


    pos[0] = R * cos(theta) * sin(phi);
    pos[1] = R * sin(theta) * sin(phi);
    pos[2] = R * cos(phi);

}


void updatePart_GF ( particle *P, gsl_rng *r, double dt, double L ) {
// differently from aGF positions are not sampled when the domain is constructed

  double deltaPos [3];

  if ( P->gf ){

    polarTransf ( deltaPos, P -> shell, gsl_rng_uniform (r), gsl_rng_uniform (r));
    //deltaPos now contains the displacements in cartesian coordinates
    P -> pos[0] += deltaPos[0];  
    P -> pos[1] += deltaPos[1];
    P -> pos[2] += deltaPos[2];
    P->shell = 0;
    P->time = P->tau_exit;

    double deltaT = trunc( (P->time+dt)/dt )*dt - P->time;

    P -> pos[0] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos[1] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos[2] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    checkBound ( P->pos, P->pos_period, L );
      P -> pos_exit[0] = P -> pos[0];
      P -> pos_exit[1] = P -> pos[1];
      P -> pos_exit[2] = P -> pos[2];
    P->time += deltaT;
    P->tau_exit = P->time;

  }
  else if ( P->burst == false ){
    // it is the case of a BM integration

    P->pos[0] = P->pos_exit[0];
    P->pos[1] = P->pos_exit[1];
    P->pos[2] = P->pos_exit[2];
    P->time = P->tau_exit;

  }

}


void updatePart_GF_P_proj ( particle *P, gsl_rng *r, double dt, double L ) {
// differently from aGF positions are not sampled when the domain is constructed

  double deltaPos [3];

  if ( P->gf  &&  P->tau_exit-P->time > (P->shell*P->shell) / P->Diff / 100 ){

    double R = drawPosNewt ( P->tau_exit-P->time, P->shell, P->Diff, gsl_rng_uniform(r) );
    double R1 = gsl_rng_uniform (r);
    double R2 = gsl_rng_uniform (r);
    polarTransf ( deltaPos, R, R1, R2 );
    //deltaPos now contains the displacements in cartesian coordinates
    P -> pos[0] += deltaPos[0];  
    P -> pos[1] += deltaPos[1];
    P -> pos[2] += deltaPos[2];      
    checkBound ( P -> pos, P -> pos_period, L );
    P->pos_exit[0] = P->pos[0];
    P->pos_exit[1] = P->pos[1];
    P->pos_exit[2] = P->pos[2];

    P->shell = 0;
    P->time = P->tau_exit;

  }
  else if (P->gf ) {

    double deltaT = P->tau_exit-P->time;
    P -> pos[0] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos[1] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos[2] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    checkBound ( P->pos, P->pos_period, L );
      P->pos_exit[0] = P->pos[0];
      P->pos_exit[1] = P->pos[1];
      P->pos_exit[2] = P->pos[2];

      P->shell = 0;
    P->time = P->tau_exit;

  }
  else if ( !P->burst ){
    // it is the case of a BM integration (no GF, no burst)

    P->pos[0] = P->pos_exit[0];
    P->pos[1] = P->pos_exit[1];
    P->pos[2] = P->pos_exit[2];
    P->time = P->tau_exit;

  }

}


void updatePart_GF_PQ_proj ( particle *P, gsl_rng *r, double dt, double L ) {
// differently from aGF, the positions are not sampled when the domain is constructed
// the function is called for the particle P that has the lowest tau_exit, while instead the tau_exitSampled can be higher

  double deltaPos[3],deltaPosFut[3];

  if ( P->gf  &&  P->tau_exit-P->time > (P->shell*P->shell) / P->Diff / 100 ){

      double R = drawPosPQ00bis ( P->tau_exit-P->time, P->tau_exitSampled-P->time, P->shell, P->Diff, gsl_rng_uniform(r) );
      double theta =  2 * M_PI * gsl_rng_uniform (r);
      double phi = acos( 2*gsl_rng_uniform (r) - 1 );

      polarTransf_angles ( deltaPos, R, theta, phi );

      P -> pos[0] += deltaPos[0];
      P -> pos[1] += deltaPos[1];
      P -> pos[2] += deltaPos[2];

      P -> totPQdispl = 1;
      P -> countPQ = 0;

      double deltaT = dt - ( P->tau_exitSampled - P->tau_exit );

      theta += gsl_ran_gaussian (r,1) * P->sqrtDiff * sqrt(2*deltaT) / R / sin(phi);
      phi += gsl_ran_gaussian (r,1) * P->sqrtDiff * sqrt(2*deltaT) / R;
      if (phi>M_PI) phi = M_PI - (phi - M_PI);
      if (phi<0) phi = abs(phi);

      polarTransf_angles ( deltaPosFut, P->shell, theta, phi );

      P->displPQ[0][0] = deltaPosFut[0] - deltaPos[0] + sqrt(2*deltaT)*P->sqrtDiff*gsl_ran_gaussian (r,1);
      P->displPQ[1][0] = deltaPosFut[1] - deltaPos[1] + sqrt(2*deltaT)*P->sqrtDiff*gsl_ran_gaussian (r,1);
      P->displPQ[2][0] = deltaPosFut[2] - deltaPos[2] + sqrt(2*deltaT)*P->sqrtDiff*gsl_ran_gaussian (r,1);

      checkBound ( P -> pos, P -> pos_period, L );
      
      P -> pos_exit [0] = P -> pos[0];
      P -> pos_exit [1] = P -> pos[1];
      P -> pos_exit [2] = P -> pos[2];

      P->shell = 0;
      P->time = P->tau_exit;
      P -> gf = false;

  }
  else if (P->gf ) {
    //this cycle is for the case of a very small time sampling (very unlikely)
    double deltaT = P->tau_exit-P->time;
    P -> pos[0] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos[1] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);
    P -> pos[2] += gsl_ran_gaussian (r,1)*P->sqrtDiff*sqrt(2*deltaT);

      checkBound ( P->pos, P->pos_period, L );
      P -> pos_exit [0] = P -> pos[0];
      P -> pos_exit [1] = P -> pos[1];
      P -> pos_exit [2] = P -> pos[2];

    P->shell = 0;
    P->time = P->tau_exit;
      P -> gf = false;

  }
  else if ( !P->burst  ){
    // it is the case of a BM integration, when the particle has been burst, it skips the update cycle

    checkBound (P->pos_exit,P->pos_period, L);

    P->pos[0] = P->pos_exit[0];
    P->pos[1] = P->pos_exit[1];
    P->pos[2] = P->pos_exit[2];
    P->time = P->tau_exit;

  }

}
