// author luigisbailo


#pragma once

void burst_P_GF ( struct particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double L ) {

  double deltaPos [3];

  // it cycles over all particles to check whether they are within the bursting radius
  for (int j=1; j<N; j++){

    int jPart = partList[j];


    if ( particles[jPart].gf==0 && distRow[j] - particles[jPart].shell < particles[iPart].burstR &&
            particles[iPart].time<particles[jPart].tau_exit){

      particles[jPart].burst = 0;
      particles[jPart].gf = 1;

      //The P function is not sampled at very small times, when the survival function S can be approximated to 1      
      if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell) /
                                                               particles[jPart].Diff/100){
      
        polarTransf ( deltaPos, drawPosNewt ( particles[iPart].time-particles[jPart].time,  particles[jPart].shell,
                                              particles[jPart].Diff, gsl_rng_uniform(r) ), gsl_rng_uniform (r), gsl_rng_uniform (r) );
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

      // "distRow[]" is updated with the new distances, and whether there is a new closest distance to insert in distRow[0] is checked
      distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L )) - particles[iPart].radius - particles[jPart].radius;
      if (distRow[j]<distRow[0]) distRow[0]=distRow[j];

    }

  }  

}


void burst_PQ_GF ( struct particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double L ) {

    double deltaPos [3];

    // it cycles over all particles to check weather they are within the bursting radius
    for (int j=1; j<N; j++){

        int jPart = partList[j];


        if ( particles[jPart].gf==0 && distRow[j] - particles[jPart].shell < particles[iPart].burstR
             && particles[iPart].time<particles[jPart].tau_exit){

            particles[jPart].burst = 0;
            particles[jPart].gf = 1;


            //The P function is not sampled at very small times, when the survival function S can be approximated to 1
            if (particles[iPart].time-particles[jPart].time> (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){

                polarTransf ( deltaPos,
                              drawPosPQ00bis ( particles[iPart].time-particles[jPart].time,
                                               particles[jPart].tau_exit-particles[jPart].time ,
                                               particles[jPart].shell, particles[jPart].Diff, gsl_rng_uniform(r) ),
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

            // "distRow[]" is updated with the new distances,
            // and weather there is a new closest distance to insert in distRow[0] is checked
            distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L ))
                          - particles[iPart].radius - particles[jPart].radius;
            if (distRow[j]<distRow[0]) distRow[0]=distRow[j];


        }

    }

}



void burst_PQ_GF_proj ( struct particle *particles, int *partList, double *distRow, gsl_rng *r,  int N, int iPart, double tau_bm, double L, double tProj ) {

    double deltaPos [3];

    // it cycles over all particles to check weather they are within the bursting radius
    for (int j=1; j<N; j++){

    int jPart = partList[j];
    //Conditions:
    // 1 - particle in GF mode
    // 2 - particle exit-time smaller than bursting time
    //     one particle may exit its domain, do the BM fractional propagation and then burst the other domain after its exit time
    // 3 - domain in within bursting distance
    if ( particles[jPart].gf == 0 &&
            particles[iPart].time<particles[jPart].tau_exitSampled &&
            distRow[j] - particles[jPart].shell < particles[iPart].burstR ){

        particles[jPart].gf = 1;

        //The P function is not sampled at very small times, when the survival function S can be approximated to 1
        if (particles[iPart].time-particles[jPart].time >
                (particles[jPart].shell*particles[jPart].shell)/particles[jPart].Diff/100){
            particles[jPart].burst = 0;

            double tempPos[3], radiusPQ;
            double theta =  2 * M_PI * gsl_rng_uniform (r);
            double phi = acos( 2*gsl_rng_uniform (r) - 1 );
            double tau_exit = particles[jPart].tau_exitSampled-particles[jPart].time;
            double t_sampling = particles[iPart].time-particles[jPart].time;
            radiusPQ = drawPosPQ00bis ( t_sampling, tau_exit, particles[jPart].shell,
                                        particles[jPart].Diff, gsl_rng_uniform(r) );

            polarTransf_angles ( deltaPos, radiusPQ, theta, phi);
            //deltaPos now contains the displacements in cartesian coordinate

            particles[jPart].pos_exit[0] = particles[jPart].pos[0] + deltaPos[0];
            particles[jPart].pos_exit[1] = particles[jPart].pos[1] + deltaPos[1];
            particles[jPart].pos_exit[2] = particles[jPart].pos[2] + deltaPos[2];
            checkBound ( particles[jPart].pos_exit, particles[jPart].pos_period, L );

            particles[jPart].countPQ = 0;

            double radius0;
            tau_exit -= t_sampling;
            radius0 = radiusPQ;
            int count_PQ = 0;


            double X_temp = deltaPos[0];
            double Y_temp = deltaPos[1];
            double Z_temp = deltaPos[2];
            double dx_trial=0;
            double dy_trial=0;
            double dz_trial=0;

            double epsilon = 0.01;

            while (tau_exit>tau_bm) {

                double Xi = gsl_rng_uniform(r);

                radiusPQ = drawPosPQbis ( tau_bm, radius0, tau_exit, particles[jPart].shell, particles[jPart].Diff, Xi);

                double R_temp = radius0;
                if (fabs(radius0-radiusPQ)<3*sqrt(6*particles[jPart].Diff*tau_bm)) {
                     do{
                        dx_trial = gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt(2 * tau_bm);
                        dy_trial = gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt(2 * tau_bm);
                        dz_trial = gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt(2 * tau_bm);

                        R_temp = sqrt(pow(X_temp + dx_trial, 2) + pow(Y_temp + dy_trial, 2) + pow(Z_temp + dz_trial, 2));

                    } while (!(R_temp > radiusPQ - epsilon &&
                             R_temp < radiusPQ + epsilon && R_temp<particles[jPart].shell) );

                }
                else{
                    if (X_temp>0){
                        theta = atan(Y_temp/X_temp);
                    }
                    else{
                        theta = atan(Y_temp/X_temp)+M_PI;

                    }
                    phi = acos (Z_temp/radius0);

                    dx_trial = (radiusPQ - radius0) * cos(theta) * sin(phi);
                    dy_trial = (radiusPQ - radius0) * sin(theta) * sin(phi);
                    dz_trial = (radiusPQ - radius0) * cos(phi);

                }
                X_temp += dx_trial;
                Y_temp += dy_trial;
                Z_temp += dz_trial;

                particles[jPart].displPQ[0][count_PQ] = dx_trial;
                particles[jPart].displPQ[1][count_PQ] = dy_trial;
                particles[jPart].displPQ[2][count_PQ] = dz_trial;

                radius0 = sqrt(X_temp * X_temp + Y_temp * Y_temp + Z_temp * Z_temp);
                tau_exit -= tau_bm;

                count_PQ ++;
                if (count_PQ>MEMORY_PQ){
                    printf("ERROR: MEMORY_PQ too small");
                    exit(EXIT_FAILURE);
                }

            }

            if (X_temp>0){
                theta = atan(Y_temp/X_temp);
            }
            else{
                theta = atan(Y_temp/X_temp)+M_PI;

            }
            phi = acos (Z_temp/radius0);

            theta += gsl_ran_gaussian (r,1) * particles[jPart].sqrtDiff * sqrt(2*tau_bm) / radius0 / sin(phi);
            phi += gsl_ran_gaussian (r,1) * particles[jPart].sqrtDiff * sqrt(2*tau_bm) / radius0;
            if (phi>M_PI) phi = M_PI - (phi - M_PI);
            if (phi<0) phi = fabs(phi);

            polarTransf_angles(tempPos, particles[jPart].shell, theta, phi);

            double deltaT = tau_bm - (particles[jPart].tau_exitSampled - particles[jPart].tau_exit);

            particles[jPart].displPQ[0][count_PQ] = tempPos[0] - X_temp +
                    sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);
            particles[jPart].displPQ[1][count_PQ] = tempPos[1] - Y_temp +
                    sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);
            particles[jPart].displPQ[2][count_PQ] = tempPos[2] - Z_temp +
                    sqrt(2*deltaT)*particles[jPart].sqrtDiff*gsl_ran_gaussian (r,1);

            particles[jPart].totPQdispl = count_PQ;

            particles[jPart].pos[0] = particles[jPart].pos_exit[0];
            particles[jPart].pos[1] = particles[jPart].pos_exit[1];
            particles[jPart].pos[2] = particles[jPart].pos_exit[2];
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
      distRow [j] = sqrt(dist2_per ( &particles[iPart], &particles[jPart], L )) -
              particles[iPart].radius - particles[jPart].radius;
      if (distRow[j]<distRow[0]) distRow[0]=distRow[j];      

    }

  }  

}