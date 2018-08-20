#pragma once

void initPos_hybGF ( particle *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm, double alpha, double L) {

    double x,y,z;

    for (int count=0; count < N_A+N_B; count++){

    if (count < N_A){

      particles[count].type = 0;
      particles[count].Diff = D_A;
      particles[count].sqrtDiff = sqrt(D_A);
      particles[count].radius = R_A;
      particles[count].R_gfrd = alpha*sqrt(D_A*tau_bm);
      particles[count].R_bd = alpha*sqrt(D_A*tau_bm);
      particles[count].burstR = alpha*sqrt(fmin(D_A,D_B)*tau_bm);

    }
    else{

      particles[count].type = 1;
      particles[count].Diff = D_B;
      particles[count].sqrtDiff = sqrt(D_B);
      particles[count].radius = R_B;
      particles[count].R_gfrd = alpha*sqrt(D_B*tau_bm);
      particles[count].R_bd = alpha*sqrt(D_A*tau_bm);
      particles[count].burstR = alpha*sqrt(fmin(D_A,D_B)*tau_bm);

    }


    particles[count].label = count;
    particles[count].time = 0;
    particles[count].tau_exit = 0;
    particles[count].tau_exitSampled=-1;
    particles[count].shell = 0; 
    particles[count].burst = false; 
    particles[count].gf = false;
    for ( int i = 0; i < MEMORY_PQ; i++){

        particles[count].displPQ[0][i] = 0;
        particles[count].displPQ[1][i] = 0;
        particles[count].displPQ[2][i] = 0;

    }
    particles[count].countPQ = 0;

    double distMin,dist;
    //The following cycle ensures that particles are not in within interaction distance at initialization
    do{

      distMin = 0;

      x = gsl_rng_uniform (r)*L;
      y = gsl_rng_uniform (r)*L;  
      z = gsl_rng_uniform (r)*L;
      
      particles[count].pos[0] = x;
      particles[count].pos[1] = y;
      particles[count].pos[2] = z;

      for ( int n=0; n<count; n++ ){

        dist = sqrt(dist2_per ( &particles[count], &particles[n], L )) - particles[count].radius - particles[n].radius; 
        if (dist<distMin) distMin = dist;

     } 

    }while (distMin<-particles[count].radius/2);

    
    particles[count].pos_exit[0] = -1;
    particles[count].pos_exit[1] = -1;
    particles[count].pos_exit[2] = -1;   
    particles[count].pos_init[0] = x;
    particles[count].pos_init[1] = y;
    particles[count].pos_init[2] = z; 
    particles[count].pos_period[0] = 0;
    particles[count].pos_period[1] = 0;
    particles[count].pos_period[2] = 0;


  }
}



void initPos_BM ( particle *particles, gsl_rng *r, int N_A, int N_B, double R_A, double R_B, double D_A, double D_B, double tau_bm, double L) {

  double x,y,z;

  for (int count=0; count < N_A+N_B; count++){

    particles[count].label = count;
    particles[count].time = 0;
    particles[count].tau_exit = 0; 
    particles[count].shell = 0; 
    particles[count].burst = false; 
    particles[count].gf = false;

    double distMin,dist;

    do{

      distMin = 0;

      x = gsl_rng_uniform (r)*L;
      y = gsl_rng_uniform (r)*L;  
      z = gsl_rng_uniform (r)*L;
      
      particles[count].pos[0] = x;
      particles[count].pos[1] = y;
      particles[count].pos[2] = z;

      for ( int n=0; n<count; n++ ){

        dist = sqrt(dist2_per ( &particles[count], &particles[n], L )) - particles[count].radius - particles[n].radius; 
        if (dist<distMin) distMin = dist;

     } 

    }while (distMin<0);

    

    particles[count].pos_exit[0] = x;
    particles[count].pos_exit[1] = y;
    particles[count].pos_exit[2] = z;   
    particles[count].pos_init[0] = x;
    particles[count].pos_init[1] = y;
    particles[count].pos_init[2] = z; 
    particles[count].pos_period[0] = 0;
    particles[count].pos_period[1] = 0;
    particles[count].pos_period[2] = 0;


    if (count < N_A){        

      particles[count].Diff = D_A;    
      particles[count].sqrtDiff = sqrt(D_A); 
      particles[count].radius = R_A;
      
    }
    else{

      particles[count].Diff = D_B;    
      particles[count].sqrtDiff = sqrt(D_B); 
      particles[count].radius = R_B;

    }
  }
}





void initShell_GF ( particle *particles, gsl_rng *r, int N, double tau_bm, double sqrt2TAU_BM, double L, int *stat ) {

  double dists [N], deltaPos[3], varPos[3];

  for (int i=0; i<N; i++){

    for (int j=0; j<N; j++) {

      if (i!=j){
        dists [j] = (sqrt (dist2_per (&particles[i],&particles[j],L) ) - particles[i].radius - particles[j].radius );
      }
      else {
        dists [j] = L;
      }

    }
    
    double R = *std::min_element( dists, dists+N ) / 2;


  if ( R > particles[i].R_bd ){

    (*stat) ++;
 
    particles[i].gf = true;
    particles[i].tau_exit += drawTimeNewt ( R, particles[i].Diff, gsl_rng_uniform(r) ); 
    particles[i].shell = R;
    particles[i].pos_exit[0]= -1;
    particles[i].pos_exit[1]= -1;
    particles[i].pos_exit[2]= -1;

  }      
  else {

    particles[i].gf = false;

    deltaPos[0] = 0;
    deltaPos[1] = 0;
    deltaPos[2] = 0;

    for ( int j=0; j<N; j++){

      if ( dists [j] < 0 && i!=j ){

        varPos[0] = particles[i].pos[0] - particles[j].pos[0];
        if (varPos[0]>L/2) varPos[0] -= L;
        else if (varPos[0]<-L/2) varPos[0] += L;
        varPos[1] = particles[i].pos[1] - particles[j].pos[1];
        if (varPos[1]>L/2) varPos[1] -= L;
        else if (varPos[1]<-L/2) varPos[1] += L;
        varPos[2] = particles[i].pos[2] - particles[j].pos[2];
        if (varPos[2]>L/2) varPos[2] -= L;
        else if (varPos[2]<-L/2) varPos[2] += L;
   
        deltaPos[0] += K*particles[i].Diff * (varPos[0]/(dists[j]+particles[i].radius+particles[j].radius)) * (-dists[j]) * tau_bm;
        deltaPos[1] += K*particles[i].Diff * (varPos[1]/(dists[j]+particles[i].radius+particles[j].radius)) * (-dists[j]) * tau_bm;
        deltaPos[2] += K*particles[i].Diff * (varPos[2]/(dists[j]+particles[i].radius+particles[j].radius)) * (-dists[j]) * tau_bm;

      }

    }
    
    particles[i].tau_exit += tau_bm;
    particles[i].pos_exit[0] = particles[i].pos[0] + deltaPos[0] + gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    particles[i].pos_exit[1] = particles[i].pos[1] + deltaPos[1] + gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    particles[i].pos_exit[2] = particles[i].pos[2] + deltaPos[2] + gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;

    checkBound (particles[i].pos_exit,particles[i].pos_period, L);

  }

 }


}


void initShell_GF_proj ( particle *particles, gsl_rng *r, int N, double tau_bm, double sqrt2TAU_BM, double L, int *stat ) {

  double dists [N], deltaPos[3], varPos[3];


  for (int i=0; i<N; i++){

    for (int j=0; j<N; j++) {

      if (i!=j){
        dists [j] = (sqrt (dist2_per (&particles[i],&particles[j],L) ) - particles[i].radius - particles[j].radius );
      }
      else {
        dists [j] = L;
      }

    }
    
   double R = *std::min_element( dists, dists+N ) / 2;


  if ( R > particles[i].R_bd ){

    (*stat) ++;
 
    particles[i].gf = true;
    particles[i].tau_exit += drawTimeNewt ( R, particles[i].Diff, gsl_rng_uniform(r) ); 
    particles[i].tau_exitSampled = particles[i].tau_exit;
    particles[i].tau_exit = trunc ( particles[i].tau_exit / tau_bm ) * tau_bm;
    particles[i].shell = R;
    particles[i].pos_exit[0]= -1;
    particles[i].pos_exit[1]= -1;
    particles[i].pos_exit[2]= -1;

  }      
  else {

    particles[i].gf = false;

    deltaPos[0] = 0;
    deltaPos[1] = 0;
    deltaPos[2] = 0;

    for ( int j=0; j<N; j++){

      if ( dists [j] < 0 && i!=j ){

        varPos[0] = particles[i].pos[0] - particles[j].pos[0];
        if (varPos[0]>L/2) varPos[0] -= L;
        else if (varPos[0]<-L/2) varPos[0] += L;
        varPos[1] = particles[i].pos[1] - particles[j].pos[1];
        if (varPos[1]>L/2) varPos[1] -= L;
        else if (varPos[1]<-L/2) varPos[1] += L;
        varPos[2] = particles[i].pos[2] - particles[j].pos[2];
        if (varPos[2]>L/2) varPos[2] -= L;
        else if (varPos[2]<-L/2) varPos[2] += L;
   
        deltaPos[0] += K*particles[i].Diff * (varPos[0]/(dists[j]+particles[i].radius+particles[j].radius)) * (-dists[j]) * tau_bm;
        deltaPos[1] += K*particles[i].Diff * (varPos[1]/(dists[j]+particles[i].radius+particles[j].radius)) * (-dists[j]) * tau_bm;
        deltaPos[2] += K*particles[i].Diff * (varPos[2]/(dists[j]+particles[i].radius+particles[j].radius)) * (-dists[j]) * tau_bm;

      }

    }
    
    particles[i].tau_exit += tau_bm;
    particles[i].pos_exit[0] = particles[i].pos[0] + deltaPos[0] + gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    particles[i].pos_exit[1] = particles[i].pos[1] + deltaPos[1] + gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
    particles[i].pos_exit[2] = particles[i].pos[2] + deltaPos[2] + gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;

    checkBound (particles[i].pos_exit,particles[i].pos_period, L);

  }

 }


}
