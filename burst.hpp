#pragma once

void burst_P_GF ( particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double L ) {

  double deltaPos [3];

  // it cycles over all particles to check weather they are within the bursting radius
  for (int j=1; j<N; j++){

    int jPart = partList[j];


    if ( particles[jPart].gf == true && distRow[j] - particles[jPart].shell < particles[iPart].burstR ){

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




void burst_PQ_GF_proj ( particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double tau_bm, double L ) {

    double deltaPos [3];

  // it cycles over all particles to check weather they are within the bursting radius
    for (int j=1; j<N; j++){

    int jPart = partList[j];

    if ( particles[jPart].gf  && distRow[j] - particles[jPart].shell < particles[iPart].burstR ){

        //        std::cout<<"Burst  " << particles[iPart].time <<std::endl;
        particles[jPart].gf = false;

        //The P function is not sampled at very small times, when the survival function S can be approximated to 1
        if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){
            particles[jPart].burst = true;

            double tempPosOld[3],tempPosNew[3], radiusPQ;
            double R1 = gsl_rng_uniform (r);
            double R2 = gsl_rng_uniform (r);
            radiusPQ = drawPosPQ00bis ( particles[iPart].time-particles[jPart].time, particles[jPart].tau_exitSampled-particles[jPart].time, particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) );

            polarTransf ( deltaPos, radiusPQ, R1, R2);
            //deltaPos now contains the displacements in cartesian coordinate

            particles[jPart].pos_exit[0] = particles[jPart].pos[0] + deltaPos[0];
            particles[jPart].pos_exit[1] = particles[jPart].pos[1] + deltaPos[1];
            particles[jPart].pos_exit[2] = particles[jPart].pos[2] + deltaPos[2];
            checkBound ( particles[jPart].pos_exit, particles[jPart].pos_period, L );


            particles[jPart].totPQdispl = int(ceil((particles[jPart].tau_exitSampled-particles[iPart].time)/tau_bm ));
            particles[jPart].countPQ = 0;

            if (particles[jPart].totPQdispl>MEMORY_PQ){
                std::cout<< "ERROR: MEMORY_PQ too small " << MEMORY_PQ << " " << particles[jPart].countPQ << std::endl;
                exit(EXIT_FAILURE);
            }

            tempPosOld[0] = deltaPos[0];
            tempPosOld[1] = deltaPos[1];
            tempPosOld[2] = deltaPos[2];

//std::cout << "Part. label: " <<  particles[jPart].label << "\t" << particles[jPart].tau_exitSampled - particles[iPart].time << "\t" << particles[jPart].totPQdispl << std::endl;

            double deltaRadius, deltaT;
            //deltaRadius is constant, it should be sampled from PQ instead
            deltaRadius = (particles[jPart].shell-radiusPQ)/(particles[jPart].totPQdispl-1);
//            std::cout << radiusPQ << "\t" << particles[jPart].shell << std::endl;
//            std::cout << particles[jPart].shell << "\t" << radiusPQ << "\t" << particles[jPart].tau_exitSampled-particles[iPart].time << "\t"  << std::endl;
            for ( int count=0; count<particles[jPart].totPQdispl-1; count++){

                //it is assumed that there is no angular displacement, i.e. always the same R1,R2 are used
                polarTransf(tempPosNew,radiusPQ+deltaRadius, R1, R2);
                particles[jPart].displPQ[0][count] = tempPosNew[0] - tempPosOld[0];
                particles[jPart].displPQ[1][count] = tempPosNew[1] - tempPosOld[1];
                particles[jPart].displPQ[2][count] = tempPosNew[2] - tempPosOld[2];
//                std::cout << count << " \t" << sqrt(pow(particles[jPart].displPQ[0][count],2) + pow(particles[jPart].displPQ[1][count],2) + pow(particles[jPart].displPQ[2][count],2)) << std::endl;
                radiusPQ += deltaRadius;
//                std::cout << radiusPQ << std::endl;
                tempPosOld[0] = tempPosNew[0];
                tempPosOld[1] = tempPosNew[1];
                tempPosOld[2] = tempPosNew[2];
            }
            if (particles[jPart].shell+0.01<radiusPQ ){
                std::cout << radiusPQ << "\t" << particles[jPart].shell << std::endl;
                std::cout<<"ERROR1: in burst PQ proj" << std::endl;
                exit(EXIT_FAILURE);
            }

            polarTransf(deltaPos, particles[jPart].shell-radiusPQ, R1, R2);
            deltaT = tau_bm;// -  (particles[jPart].tau_exitSampled - particles[jPart].tau_exit);

            particles[jPart].displPQ[0][particles[jPart].totPQdispl-1] = sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1); // + random displacement after exiting the old shell ;
            particles[jPart].displPQ[1][particles[jPart].totPQdispl-1] = sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);
            particles[jPart].displPQ[2][particles[jPart].totPQdispl-1] = sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);


//            particles[jPart].displPQ[0][particles[jPart].totPQdispl-1] = deltaPos[0] - tempPosNew[0] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1); // + random displacement after exiting the old shell ;
//            particles[jPart].displPQ[1][particles[jPart].totPQdispl-1] = deltaPos[1] - tempPosNew[1] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);
//            particles[jPart].displPQ[2][particles[jPart].totPQdispl-1] = deltaPos[2] - tempPosNew[2] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);
//            std::cout << particles[jPart].totPQdispl-1 << " \t" << sqrt(pow(particles[jPart].displPQ[0][particles[jPart].totPQdispl-1],2) + pow(particles[jPart].displPQ[1][particles[jPart].totPQdispl-1],2) + pow(particles[jPart].displPQ[2][particles[jPart].totPQdispl-1],2)) << std::endl;

            particles[jPart].pos[0] = particles[jPart].pos_exit[0];
            particles[jPart].pos[1] = particles[jPart].pos_exit[1];
            particles[jPart].pos[2] = particles[jPart].pos_exit[2];
            particles[jPart].shell = 0;
            particles[jPart].time = particles[iPart].time;
            particles[jPart].tau_exit = particles[iPart].time;

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


