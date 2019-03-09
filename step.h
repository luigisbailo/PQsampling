// author luigisbailo


#pragma once

void GFstep_GF ( struct particle *myPart, gsl_rng *r, double R ){
//The shell radius "R" has been already determined with getR () and it has been already checked that the domain can be constructed
    //extraction of the exit time and exit position in polar coordinates
    myPart->tau_exit += drawTimeNewt ( R, myPart->Diff, gsl_rng_uniform(r) );
    myPart->tau_exitSampled = myPart->tau_exit;
    myPart->shell = R;

}


void GFstep_GF_proj ( struct particle *myPart, gsl_rng *r, double R, double dt){
//The shell radius "R" has been already determined with getR () and it has been already checked that the domain can be constructed

  //extraction of the exit time and exit position in polar coordinates
  myPart->tau_exit += drawTimeNewt ( R, myPart->Diff, gsl_rng_uniform(r) );
  myPart->tau_exitSampled = myPart -> tau_exit;
  myPart->tau_exit = trunc( myPart->tau_exit / dt ) * dt;
  myPart->shell = R;

    
}


void BMstep ( struct particle *particles, int *partList, double *distRow, gsl_rng *r, double tau_bm,
              double sqrt2TAU_BM, int N, double L ) {

  double dist,deltaPos[3], varPos[3];

  deltaPos[0] = 0;
  deltaPos[1] = 0;
  deltaPos[2] = 0;

  for( int j =1; j<N; j ++){

    int jPart = partList[j];

    if (particles[jPart].gf==0) continue;

    dist = distRow[j];

    if ( dist < 0 ){
      // if the distance to another particle is lower than R_INTER we take into account their interaction

      varPos[0] = particles[partList[0]].pos[0] - particles[jPart].pos[0];
        if (varPos[0]>L/2) varPos[0] -= L;
        else if (varPos[0]<-L/2) varPos[0] += L;
      varPos[1] = particles[partList[0]].pos[1] - particles[jPart].pos[1];
        if (varPos[1]>L/2) varPos[1] -= L;
        else if (varPos[1]<-L/2) varPos[1] += L;
      varPos[2] = particles[partList[0]].pos[2] - particles[jPart].pos[2];
        if (varPos[2]>L/2) varPos[2] -= L;
        else if (varPos[2]<-L/2) varPos[2] += L;

      // The particles interact via a repulsive harmonic interaction
      deltaPos[0] += K*particles[partList[0]].Diff *
              (varPos[0]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;
      deltaPos[1] += K*particles[partList[0]].Diff *
              (varPos[1]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;
      deltaPos[2] += K*particles[partList[0]].Diff *
              (varPos[2]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;

    }

  }

  particles[partList[0]].pos_exit[0] +=  gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff *
                                                 sqrt2TAU_BM + deltaPos[0];
  particles[partList[0]].pos_exit[1] +=  gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff *
                                                 sqrt2TAU_BM + deltaPos[1];
  particles[partList[0]].pos_exit[2] +=  gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff *
                                                 sqrt2TAU_BM + deltaPos[2];


  checkBound (particles[partList[0]].pos_exit, particles[partList[0]].pos_period, L );

    particles[partList[0]].tau_exit += tau_bm;
}


void BMstepPQ ( struct particle *particles, int *partList, double *distRow, gsl_rng *r,
                double tau_bm, double sqrt2TAU_BM, int N, double L ) {

    double dist,deltaPosInt[3], deltaPosDiff[3], varPos[3];

    //diffusive displacement
    deltaPosDiff[0] = 0;
    deltaPosDiff[1] = 0;
    deltaPosDiff[2] = 0;

    //interaction displacement
    deltaPosInt[0] = 0;
    deltaPosInt[1] = 0;
    deltaPosInt[2] = 0;

    for( int j=1; j<N; j ++){

        int jPart = partList[j];

        if (particles[jPart].gf) continue;

        dist = distRow[j];

        if ( dist < 0 ){
            // if the distance with another particle is lower than R_INTER we take into account their interaction

            varPos[0] = particles[partList[0]].pos[0] - particles[jPart].pos[0];
            if (varPos[0]>L/2) varPos[0] -= L;
            else if (varPos[0]<-L/2) varPos[0] += L;
            varPos[1] = particles[partList[0]].pos[1] - particles[jPart].pos[1];
            if (varPos[1]>L/2) varPos[1] -= L;
            else if (varPos[1]<-L/2) varPos[1] += L;
            varPos[2] = particles[partList[0]].pos[2] - particles[jPart].pos[2];
            if (varPos[2]>L/2) varPos[2] -= L;
            else if (varPos[2]<-L/2) varPos[2] += L;

            // The particles interact via a repulsive harmonic interaction
            deltaPosInt[0] += K*particles[partList[0]].Diff *
                    (varPos[0]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;
            deltaPosInt[1] += K*particles[partList[0]].Diff *
                    (varPos[1]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;
            deltaPosInt[2] += K*particles[partList[0]].Diff *
                    (varPos[2]/(dist+particles[partList[0]].radius+particles[jPart].radius)) * (-dist) * tau_bm;

        }

    }

    if ( particles[partList[0]].tau_exitSampled>particles[partList[0]].time) {

        deltaPosDiff[0] = particles[partList[0]].displPQ[0][particles[partList[0]].countPQ];
        deltaPosDiff[1] = particles[partList[0]].displPQ[1][particles[partList[0]].countPQ];
        deltaPosDiff[2] = particles[partList[0]].displPQ[2][particles[partList[0]].countPQ];
        particles[partList[0]].countPQ ++;
    }
    else{

        deltaPosDiff[0] = gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt2TAU_BM;
        deltaPosDiff[1] = gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt2TAU_BM;
        deltaPosDiff[2] = gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt2TAU_BM;
     }

    particles[partList[0]].pos_exit[0] = particles[partList[0]].pos[0] + deltaPosInt[0] + deltaPosDiff[0];
    particles[partList[0]].pos_exit[1] = particles[partList[0]].pos[1] + deltaPosInt[1] + deltaPosDiff[1];
    particles[partList[0]].pos_exit[2] = particles[partList[0]].pos[2] + deltaPosInt[2] + deltaPosDiff[2];
    particles[partList[0]].tau_exit += tau_bm;

}


void synchPart_P_GF ( struct particle *particles, int *partList, gsl_rng *r, int N, double Tsynch, double L ) {


  double synchPos [3];

  for ( int n=0; n<N; n++) {

    if (particles[n].active == 1)
        continue;

    if ( Tsynch < particles[n].time && particles[n].active == 0 ){

       printf("ERROR: synch\n");
       exit (EXIT_FAILURE);

    }


    particles[n].gf = 0;
    particles[n].burst = 0;

    if ( particles[n].shell>0  && Tsynch-particles[n].time> (particles[n].shell*particles[n].shell)/particles[n].Diff/100 ){

      double Rsynch =  drawPosNewt ( Tsynch-particles[n].time, particles[n].shell, particles[n].Diff, gsl_rng_uniform(r) );
      polarTransf ( synchPos, Rsynch, gsl_rng_uniform (r), gsl_rng_uniform (r));

      particles[n].pos[0] += synchPos[0];
      particles[n].pos[1] += synchPos[1];
      particles[n].pos[2] += synchPos[2];
      checkBound (particles[n].pos,particles[n].pos_period, L );
      particles[n].pos_exit[0] = particles[n].pos[0];
      particles[n].pos_exit[1] = particles[n].pos[1];
      particles[n].pos_exit[2] = particles[n].pos[2];
      particles[n].shell = 0;
      particles[n].time = Tsynch;
      particles[n].tau_exit = Tsynch;
      particles[n].tau_exitSampled = Tsynch;

    }
    else if (particles[n].shell>0) {

      particles[n].pos[0] += gsl_ran_gaussian (r,1) * particles[n].sqrtDiff * sqrt( 2 * ( Tsynch-particles[n].time ) ) ;
      particles[n].pos[1] += gsl_ran_gaussian (r,1) * particles[n].sqrtDiff * sqrt( 2 * ( Tsynch-particles[n].time ) ) ;
      particles[n].pos[2] += gsl_ran_gaussian (r,1) * particles[n].sqrtDiff * sqrt( 2 * ( Tsynch-particles[n].time ) ) ;
      checkBound (particles[n].pos,particles[n].pos_period, L );
      particles[n].pos_exit[0] = particles[n].pos[0];
      particles[n].pos_exit[1] = particles[n].pos[1];
      particles[n].pos_exit[2] = particles[n].pos[2];
      particles[n].shell = 0;
      particles[n].time = Tsynch;
      particles[n].tau_exit = Tsynch;
      particles[n].tau_exitSampled = Tsynch;


    }
    else{
      particles[n].pos[0] = particles[n].pos_exit[0];
      particles[n].pos[1] = particles[n].pos_exit[1];
      particles[n].pos[2] = particles[n].pos_exit[2];
      particles[n].shell = 0;
      particles[n].time = Tsynch;
      particles[n].tau_exit = Tsynch;
      particles[n].tau_exitSampled = Tsynch;

    }

  }

}


void synchPart_PQ_GF ( struct particle *particles, int *partList, gsl_rng *r, int N, double Tsynch, double L ) {


  double synchPos [3];

  for ( int n=0; n<N; n++) {

      if (!particles[n].active)
          continue;

    particles[n].gf = 0;
    particles[n].burst = 0;
    if (particles[n].time>=Tsynch){
        continue;
    }
    if (particles[n].shell>0 &&
            Tsynch-particles[n].time > (particles[n].shell*particles[n].shell)/particles[n].Diff/100  ){

        double Rsynch =  drawPosPQ00bis( Tsynch-particles[n].time, particles[n].tau_exitSampled-particles[n].time,
                                         particles[n].shell, particles[n].Diff, gsl_rng_uniform(r) );
        polarTransf ( synchPos, Rsynch, gsl_rng_uniform (r), gsl_rng_uniform (r));
        particles[n].pos[0] += synchPos[0];
        particles[n].pos[1] += synchPos[1];
        particles[n].pos[2] += synchPos[2];

    }
    else if (particles[n].tau_exitSampled > particles[n].time){

        particles[n].pos[0] += particles[n].displPQ[0][particles[n].countPQ];
        particles[n].pos[1] += particles[n].displPQ[1][particles[n].countPQ];
        particles[n].pos[2] += particles[n].displPQ[2][particles[n].countPQ];
        particles[n].countPQ ++;

    }
    else{
        //It is the case of BM
        particles[n].pos[0] += gsl_ran_gaussian (r,1) * particles[n].sqrtDiff * sqrt( 2 * ( Tsynch-particles[n].time ) );
        particles[n].pos[1] += gsl_ran_gaussian (r,1) * particles[n].sqrtDiff * sqrt( 2 * ( Tsynch-particles[n].time ) );
        particles[n].pos[2] += gsl_ran_gaussian (r,1) * particles[n].sqrtDiff * sqrt( 2 * ( Tsynch-particles[n].time ) );

      }
      checkBound (particles[n].pos, particles[n].pos_period, L );

      particles[n].pos_exit[0] = particles[n].pos[0];
      particles[n].pos_exit[1] = particles[n].pos[1];
      particles[n].pos_exit[2] = particles[n].pos[2];

      particles[n].shell = 0;
      particles[n].time = Tsynch;
      particles[n].tau_exit = Tsynch;

  }
    
}