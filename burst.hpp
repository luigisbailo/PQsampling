#pragma once

void burst_P_GF ( particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double L ) {

  double deltaPos [3];

  // it cycles over all particles to check weather they are within the bursting radius
  for (int j=1; j<N; j++){

    int jPart = partList[j];


    if ( particles[jPart].gf && distRow[j] - particles[jPart].shell < particles[iPart].burstR && particles[iPart].time<particles[jPart].tau_exit){
//
      particles[jPart].burst = true;  
      particles[jPart].gf = false;

//        std::cout << particles[jPart].label << std::endl ;
//        		 std::cout << std::setprecision(6);
//		 printPos_per ( particles, partList, N );
//		 // printDist_per (particles, partList, N, L);
//		 std::cout << "\n";


      //The P function is not sampled at very small times, when the survival function S can be approximated to 1      
      if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){
      
        polarTransf ( deltaPos, drawPosNewt ( particles[iPart].time-particles[jPart].time,  particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) ),
                      gsl_rng_uniform (r), gsl_rng_uniform (r) );
        //deltaPos now contains the displacements in cartesian coordinates
        
        particles[jPart].pos[0] += deltaPos[0];  
        particles[jPart].pos[1] += deltaPos[1];
        particles[jPart].pos[2] += deltaPos[2];      
        checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
        particles[jPart].pos_exit[0] = particles[jPart].pos[0];
        particles[jPart].pos_exit[1] = particles[jPart].pos[1];
        particles[jPart].pos_exit[2] = particles[jPart].pos[2];
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
          particles[jPart].pos_exit[0] = particles[jPart].pos[0];
          particles[jPart].pos_exit[1] = particles[jPart].pos[1];
          particles[jPart].pos_exit[2] = particles[jPart].pos[2];
        particles[jPart].shell = 0;
        particles[jPart].time = particles[iPart].time;
        particles[jPart].tau_exit = particles[iPart].time;

      }
      else{
        //In case the domain is burst at the same time of the construction

        particles[jPart].shell = 0;
        particles[jPart].tau_exit = particles[iPart].time;

      }

      // "distRow[]" is updated with the new distances, and weather there is a new closest distance to insert in distRow[0] is checked    
      distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L )) - particles[iPart].radius - particles[jPart].radius;
      if (distRow[j]<distRow[0]) distRow[0]=distRow[j];

//        std::cout << std::setprecision(6);
//        printPos_per ( particles, partList, N );
//        // printDist_per (particles, partList, N, L);
//        std::cout << "\n";

    }

  }  

}


void burst_PQ_GF ( particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double L ) {

    double deltaPos [3];

    // it cycles over all particles to check weather they are within the bursting radius
    for (int j=1; j<N; j++){

        int jPart = partList[j];


        if ( particles[jPart].gf && distRow[j] - particles[jPart].shell < particles[iPart].burstR && particles[iPart].time<particles[jPart].tau_exit){
//
            particles[jPart].burst = true;
            particles[jPart].gf = false;

//        std::cout << particles[jPart].label << std::endl ;
//        		 std::cout << std::setprecision(6);
//		 printPos_per ( particles, partList, N );
//		 // printDist_per (particles, partList, N, L);
//		 std::cout << "\n";


            //The P function is not sampled at very small times, when the survival function S can be approximated to 1
            if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){

                polarTransf ( deltaPos, drawPosPQ00bis ( particles[iPart].time-particles[jPart].time, particles[jPart].tau_exit-particles[jPart].time ,  particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) ),
                              gsl_rng_uniform (r), gsl_rng_uniform (r) );
                //deltaPos now contains the displacements in cartesian coordinates

                particles[jPart].pos[0] += deltaPos[0];
                particles[jPart].pos[1] += deltaPos[1];
                particles[jPart].pos[2] += deltaPos[2];
                checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
                particles[jPart].pos_exit[0] = particles[jPart].pos[0];
                particles[jPart].pos_exit[1] = particles[jPart].pos[1];
                particles[jPart].pos_exit[2] = particles[jPart].pos[2];
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
                particles[jPart].pos_exit[0] = particles[jPart].pos[0];
                particles[jPart].pos_exit[1] = particles[jPart].pos[1];
                particles[jPart].pos_exit[2] = particles[jPart].pos[2];
                particles[jPart].shell = 0;
                particles[jPart].time = particles[iPart].time;
                particles[jPart].tau_exit = particles[iPart].time;

            }
            else{
                //In case the domain is burst at the same time of the construction

                particles[jPart].shell = 0;
                particles[jPart].tau_exit = particles[iPart].time;

            }

            // "distRow[]" is updated with the new distances, and weather there is a new closest distance to insert in distRow[0] is checked
            distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L )) - particles[iPart].radius - particles[jPart].radius;
            if (distRow[j]<distRow[0]) distRow[0]=distRow[j];

//        std::cout << std::setprecision(6);
//        printPos_per ( particles, partList, N );
//        // printDist_per (particles, partList, N, L);
//        std::cout << "\n";

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
            double theta =  2 * M_PI * gsl_rng_uniform (r);
            double phi = acos( 2*gsl_rng_uniform (r) - 1 );
            double tau_exit = particles[jPart].tau_exitSampled-particles[jPart].time;
            double t_sampling = particles[iPart].time-particles[jPart].time;
            radiusPQ = drawPosPQ00bis ( t_sampling, tau_exit, particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) );

            polarTransf_angles ( deltaPos, radiusPQ, theta, phi);
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

            double deltaXtot = deltaPos[0];
            double deltaYtot = deltaPos[1];
            double deltaZtot = deltaPos[2];

            double corrCoeff = 1/sqrt(particles[jPart].Diff * tau_bm);

            while (tau_exit>tau_bm) {

                polarTransf_angles(tempPosOld, radius0, theta, phi);

                radiusPQ = drawPosPQbis ( tau_bm, radius0, tau_exit, particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r));
//                double deltaR = radiusPQ - radius0;
//
//                double dx = gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt(2*tau_bm);
//                double dy = gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt(2*tau_bm);


                theta += gsl_ran_gaussian (r,1) * particles[jPart].sqrtDiff * sqrt(2*tau_bm) / radius0 / sin(phi)/(1+1/radius0/corrCoeff);
                phi += gsl_ran_gaussian (r,1) * particles[jPart].sqrtDiff * sqrt(2*tau_bm) / radius0/(1+1/radius0/corrCoeff);
                if (phi>M_PI) phi = M_PI - (phi - M_PI);
                if (phi<0) phi = abs(phi);
//                std::cout << theta << "\t" << phi << std::endl;


                polarTransf_angles(tempPosNew, radiusPQ, theta, phi);
//
//                deltaXtot += tempPosNew [0];
//                deltaYtot += tempPosNew [1];
//                deltaZtot += tempPosNew [2];
//
//                radius0 = sqrt(pow(deltaXtot,2) + pow(deltaYtot,2) + pow(deltaZtot,2));


//                double deltax =  gsl_ran_gaussian(r, 1) * particles[jPart].sqrtDiff * sqrt(2*tau_bm);
//                double deltay =  gsl_ran_gaussian(r, 1) * particles[jPart].sqrtDiff * sqrt(2*tau_bm);



//                if (theta>2*M_PI) theta -= 2*M_PI;
//                if (theta<0) theta = 2*M_PI + theta;



                //it is assumed that there is no angular displacement, i.e. always the same R1,R2 are used
                particles[jPart].displPQ[0][count_PQ] = tempPosNew[0] - tempPosOld[0];
                particles[jPart].displPQ[1][count_PQ] = tempPosNew[1] - tempPosOld[1];
                particles[jPart].displPQ[2][count_PQ] = tempPosNew[2] - tempPosOld[2];
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
//std::cout << std::endl;

            polarTransf_angles(tempPosOld, radius0, theta, phi);

            theta += gsl_ran_gaussian (r,1) * particles[jPart].sqrtDiff * sqrt(2*tau_bm) / radius0 / sin(phi)/(1+1/radius0/corrCoeff);   //20000000/radius0/radius0;
            phi += gsl_ran_gaussian (r,1) * particles[jPart].sqrtDiff * sqrt(2*tau_bm) / radius0/(1+1/radius0/corrCoeff); //20000000/radius0/radius0;
            if (phi>M_PI) phi = M_PI - (phi - M_PI);
            if (phi<0) phi = abs(phi);

            polarTransf_angles(tempPosNew, particles[jPart].shell, theta, phi);

            double deltaT = tau_bm -  (particles[jPart].tau_exitSampled - particles[jPart].tau_exit);

            particles[jPart].displPQ[0][count_PQ] = tempPosNew[0] - tempPosOld[0] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1); // + random displacement after exiting the old shell ;
            particles[jPart].displPQ[1][count_PQ] = tempPosNew[1] - tempPosOld[1] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);
            particles[jPart].displPQ[2][count_PQ] = tempPosNew[2] - tempPosOld[2] + sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);

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
//
//void burst_PQ_GF ( particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double L ) {
//
//    double deltaPos [3];
//
//    // it cycles over all particles to check weather they are within the bursting radius
//    for (int j=1; j<N; j++){
//
//        int jPart = partList[j];
//
//
//        if ( particles[jPart].gf  && distRow[j] - particles[jPart].shell < particles[iPart].burstR ){
////            printPos_per ( particles, partList, N );
//
//            particles[jPart].burst = true;
//            particles[jPart].gf = false;
//
//            if (particles[iPart].time>particles[jPart].tau_exit){
//                //This can happen for fractional propagation after the domain exit
//                particles[jPart].burst = false;
//                continue;
//            }
//            //The P function is not sampled at very small times, when the survival function S can be approximated to 1
//            else if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){
////std::cout << particles[jPart].tau_exit-particles[jPart].time  << "\t" << particles[iPart].time-particles[jPart].time << std::endl;
//
//                polarTransf ( deltaPos, drawPosPQ00bis ( particles[iPart].time-particles[jPart].time, particles[jPart].tau_exit-particles[jPart].time ,  particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) ),
//                              gsl_rng_uniform (r), gsl_rng_uniform (r) );
//                //deltaPos now contains the displacements in cartesian coordinates
////std::cout << deltaPos[0] << "\t" << deltaPos[1] << "\t" << deltaPos[2] << std::endl;
//                particles[jPart].pos[0] += deltaPos[0];
//                particles[jPart].pos[1] += deltaPos[1];
//                particles[jPart].pos[2] += deltaPos[2];
//                checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
//                particles[jPart].pos_exit[0] = particles[jPart].pos[0];
//                particles[jPart].pos_exit[1] = particles[jPart].pos[1];
//                particles[jPart].pos_exit[2] = particles[jPart].pos[2];
//                particles[jPart].shell = 0;
//                particles[jPart].time = particles[iPart].time;
//                particles[jPart].tau_exit = particles[iPart].time;
//
//            }
//            else if (particles[iPart].time>particles[jPart].time){
//
//                //At very small times, the bursting procedure consists simply in a brownian motion integration step
//                double sqrt2dt = sqrt (2*(particles[iPart].time-particles[jPart].time));
//                particles[jPart].pos[0] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
//                particles[jPart].pos[1] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
//                particles[jPart].pos[2] += gsl_ran_gaussian (r,1)*particles[jPart].sqrtDiff * sqrt2dt;
//                checkBound ( particles[jPart].pos, particles[jPart].pos_period, L );
//                particles[jPart].pos_exit[0] = particles[jPart].pos[0];
//                particles[jPart].pos_exit[1] = particles[jPart].pos[1];
//                particles[jPart].pos_exit[2] = particles[jPart].pos[2];
//                particles[jPart].shell = 0;
//                particles[jPart].time = particles[iPart].time;
//                particles[jPart].tau_exit = particles[iPart].time;
//
//            }
//            else{
//
//                //In case the domain is burst at the same time of the construction
//                 particles[jPart].shell = 0;
//                 particles[jPart].tau_exit = particles[iPart].time;
//
//            }
//
//            // "distRow[]" is updated with the new distances, and weather there is a new closest distance to insert in distRow[0] is checked
//            distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L )) - particles[iPart].radius - particles[jPart].radius;
//            if (distRow[j]<distRow[0]) distRow[0]=distRow[j];
//
//
//        }
//
//    }
//
//}
//
//
