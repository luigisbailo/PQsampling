#pragma once

void burst_P_GF ( particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double L ) {

  double deltaPos [3];

  // it cycles over all particles to check weather they are within the bursting radius
  for (int j=1; j<N; j++){

    int jPart = partList[j];


    if ( particles[jPart].gf && distRow[j] - particles[jPart].shell < particles[iPart].burstR ){

      particles[jPart].burst = true;  
      particles[jPart].gf = false;

//        std::cout<<"Burst  " << particles[iPart].time <<std::endl;

      //The P function is not sampled at very small times, when the survival function S can be approximated to 1      
      if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){
      
        polarTransf ( deltaPos, drawPosNewt ( particles[iPart].time-particles[jPart].time,  particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) ),
                      gsl_rng_uniform (r), gsl_rng_uniform (r) );
        //deltaPos now contains the displacements in cartesian coordinates
        
        particles[jPart].pos[0] += deltaPos[0];  
        particles[jPart].pos[1] += deltaPos[1];
        particles[jPart].pos[2] += deltaPos[2];      
        checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
        particles[jPart].pos_exit[0] = -1;
        particles[jPart].pos_exit[1] = -1;
        particles[jPart].pos_exit[2] = -1;
        particles[jPart].shell = 0;
        particles[jPart].time = particles[iPart].time;
        particles[jPart].tau_exit = particles[iPart].time;

      }
      else if (particles[iPart].time>particles[jPart].time){

        //At very small times, the bursting procedure consists simply in a brownian motion integration step 
        double sqrt2dt = sqrt (2*(particles[iPart].time-particles[jPart].time));
        particles[jPart].pos[0] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        particles[jPart].pos[1] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        particles[jPart].pos[2] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
        particles[jPart].shell = 0;
        particles[jPart].time = particles[iPart].time;
        particles[jPart].tau_exit = particles[iPart].time;

      }
      else{
        
        //In case the domain is burst at the same time of the construction
        particles[jPart].pos_exit[0]=-1;
        particles[jPart].pos_exit[1]=-1;
        particles[jPart].pos_exit[2]=-1;
        particles[jPart].shell = 0;
        particles[jPart].tau_exit = particles[iPart].time;

      }

      // "distRow[]" is updated with the new distances, and weather there is a new closest distance to insert in distRow[0] is checked    
      distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L )) - particles[iPart].radius - particles[jPart].radius;
      if (distRow[j]<distRow[0]) distRow[0]=distRow[j];      
 
    
    }

  }  

}




void burst_PQ_GF_proj ( particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double tau_bm, double L, double tProj ) {

    double deltaPos [3];

  // it cycles over all particles to check weather they are within the bursting radius
    for (int j=1; j<N; j++){

    int jPart = partList[j];

    if ( particles[jPart].gf  && distRow[j] - particles[jPart].shell < particles[iPart].burstR ){

//                std::cout<<"Burst  " << particles[iPart].time <<std::endl;
        particles[jPart].gf = false;

        //The P function is not sampled at very small times, when the survival function S can be approximated to 1
        if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){
            particles[jPart].burst = true;

            double tempPosOld[3],tempPosNew[3], radiusPQ;
            double R1 = gsl_rng_uniform (r);
            double R2 = gsl_rng_uniform (r);
            double tau_exit = particles[jPart].tau_exitSampled-particles[jPart].time;
            double t_sampling = particles[iPart].time-particles[jPart].time;
            radiusPQ = drawPosPQ00bis ( t_sampling, tau_exit, particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) );

            polarTransf ( deltaPos, radiusPQ, R1, R2);
            //deltaPos now contains the displacements in cartesian coordinate

            particles[jPart].pos_exit[0] = particles[jPart].pos[0] + deltaPos[0];
            particles[jPart].pos_exit[1] = particles[jPart].pos[1] + deltaPos[1];
            particles[jPart].pos_exit[2] = particles[jPart].pos[2] + deltaPos[2];
            checkBound ( particles[jPart].pos_exit, particles[jPart].pos_period, L );

            particles[jPart].countPQ = 0;

            tempPosOld[0] = deltaPos[0];
            tempPosOld[1] = deltaPos[1];
            tempPosOld[2] = deltaPos[2];

            double radius0;
            tau_exit -= t_sampling;
            radius0 = radiusPQ;
            int count_PQ = 0;

            double totDispl[3];
            totDispl[0]=0;
            totDispl[1]=0;
            totDispl[2]=0;
            while (tau_exit>tau_bm) {
//
//                R1 += gsl_ran_gaussian (r,1)/1000/radius0;
//                R2 += gsl_ran_gaussian (r,1)/1000/radius0;
//                if (R1<0) R1 = abs(R1);
//                if (R2<0) R2 = abs(R2);
                radiusPQ = drawPosPQbis ( tau_bm, radius0, tau_exit, particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r));
                polarTransf(tempPosNew, radiusPQ-radius0, R1, R2);
                //it is assumed that there is no angular displacement, i.e. always the same R1,R2 are used
                particles[jPart].displPQ[0][count_PQ] = tempPosNew[0];
                particles[jPart].displPQ[1][count_PQ] = tempPosNew[1];
                particles[jPart].displPQ[2][count_PQ] = tempPosNew[2];
                tempPosOld[0] = tempPosNew[0];
                tempPosOld[1] = tempPosNew[1];
                tempPosOld[2] = tempPosNew[2];
                tau_exit -= tau_bm;
                radius0 = radiusPQ;

                count_PQ ++;
                if (count_PQ>MEMORY_PQ){
                    std::cout<< "ERROR: MEMORY_PQ too small " << MEMORY_PQ << " " << particles[jPart].countPQ << std::endl;
                    exit(EXIT_FAILURE);
                }
                totDispl[0] += particles[jPart].displPQ[0][count_PQ];
                totDispl[1] += particles[jPart].displPQ[1][count_PQ];
                totDispl[2] += particles[jPart].displPQ[2][count_PQ];

            }
            polarTransf(tempPosNew, particles[jPart].shell-radius0, R1, R2);

            double deltaT = tau_bm -  (particles[jPart].tau_exitSampled - particles[jPart].tau_exit);

            particles[jPart].displPQ[0][count_PQ] = tempPosNew[0] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1); // + random displacement after exiting the old shell ;
            particles[jPart].displPQ[1][count_PQ] = tempPosNew[1] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);
            particles[jPart].displPQ[2][count_PQ] = tempPosNew[2] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);

            particles[jPart].totPQdispl = count_PQ;

            particles[jPart].pos[0] = particles[jPart].pos_exit[0];
            particles[jPart].pos[1] = particles[jPart].pos_exit[1];
            particles[jPart].pos[2] = particles[jPart].pos_exit[2];
            particles[jPart].shell = 0;
            particles[jPart].time = particles[iPart].time;
            particles[jPart].tau_exit = particles[iPart].time;
            if (particles[jPart].tau_exitSampled>tProj)
                particles[jPart].tau_exitSampled = particles[iPart].time - tau_bm;
            // tau_exitSampled remains unchanged

        }
      else if (particles[iPart].time>particles[jPart].time){

        //At very small times, the bursting procedure consists simply in a brownian motion integration step 
        double sqrt2dt = sqrt (2*(particles[iPart].time-particles[jPart].time));
        particles[jPart].pos[0] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        particles[jPart].pos[1] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        particles[jPart].pos[2] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
        checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
        particles[jPart].pos_exit[0]=particles[jPart].pos[0];
        particles[jPart].pos_exit[1]=particles[jPart].pos[1];
        particles[jPart].pos_exit[2]=particles[jPart].pos[2];
        particles[jPart].shell = 0;
        particles[jPart].time = particles[iPart].time;
        particles[jPart].tau_exit = particles[iPart].time;
        particles[jPart].tau_exitSampled = particles[iPart].time;

      }
      else{
        
        //In case the domain is burst at the same time of the construction
        particles[jPart].pos_exit[0]=particles[jPart].pos[0];
        particles[jPart].pos_exit[1]=particles[jPart].pos[1];
        particles[jPart].pos_exit[2]=particles[jPart].pos[2];
        particles[jPart].shell = 0;
        particles[jPart].tau_exit = particles[iPart].time;
        particles[jPart].tau_exitSampled = particles[iPart].time;

        }

      // "distRow[]" is updated with the new distances, and weather there is a new closest distance to insert in distRow[0] is checked    
      distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L )) - particles[iPart].radius - particles[jPart].radius;
      if (distRow[j]<distRow[0]) distRow[0]=distRow[j];      

    }

  }  

}

void burst_PQ_GF ( particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double L ) {

    double deltaPos [3];

    // it cycles over all particles to check weather they are within the bursting radius
    for (int j=1; j<N; j++){

        int jPart = partList[j];


        if ( particles[jPart].gf  && distRow[j] - particles[jPart].shell < particles[iPart].burstR ){

            particles[jPart].burst = true;
            particles[jPart].gf = false;

            //The P function is not sampled at very small times, when the survival function S can be approximated to 1
            if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){

                polarTransf ( deltaPos, drawPosPQ00bis ( particles[iPart].time-particles[jPart].time, particles[jPart].tau_exit-particles[jPart].time ,  particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) ),
                              gsl_rng_uniform (r), gsl_rng_uniform (r) );
                //deltaPos now contains the displacements in cartesian coordinates

                particles[jPart].pos[0] += deltaPos[0];
                particles[jPart].pos[1] += deltaPos[1];
                particles[jPart].pos[2] += deltaPos[2];
                checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
                particles[jPart].pos_exit[0] = -1;
                particles[jPart].pos_exit[1] = -1;
                particles[jPart].pos_exit[2] = -1;
                particles[jPart].shell = 0;
                particles[jPart].time = particles[iPart].time;
                particles[jPart].tau_exit = particles[iPart].time;

            }
            else if (particles[iPart].time>particles[jPart].time){

                //At very small times, the bursting procedure consists simply in a brownian motion integration step
                double sqrt2dt = sqrt (2*(particles[iPart].time-particles[jPart].time));
                particles[jPart].pos[0] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
                particles[jPart].pos[1] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
                particles[jPart].pos[2] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
                checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
                particles[jPart].shell = 0;
                particles[jPart].time = particles[iPart].time;
                particles[jPart].tau_exit = particles[iPart].time;

            }
            else{

                //In case the domain is burst at the same time of the construction
                particles[jPart].pos_exit[0]=-1;
                particles[jPart].pos_exit[1]=-1;
                particles[jPart].pos_exit[2]=-1;
                particles[jPart].shell = 0;
                particles[jPart].tau_exit = particles[iPart].time;

            }

            // "distRow[]" is updated with the new distances, and weather there is a new closest distance to insert in distRow[0] is checked
            distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L )) - particles[iPart].radius - particles[jPart].radius;
            if (distRow[j]<distRow[0]) distRow[0]=distRow[j];


        }

    }

}


