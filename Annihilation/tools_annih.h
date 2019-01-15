void BMstep_annih ( particle *particles, int *partList, double *distRow, gsl_rng *r, double tau_bm, double sqrt2TAU_BM, int N, double L, double Tsim ) {

//sqrtTAU_BM is sqrt(2*TAU_BM)

    double dist;

    for( int j =1; j<N; j ++){

        int jPart = partList[j];

        if (particles[jPart].gf | !particles[jPart].active) continue;

        dist = distRow[j];

        if ( dist < 0 && particles[partList[0]].type==particles[jPart].type ){
            // if the distance with another particle is lower than R_INTER we take into account their interaction
            // varPos is the cartesian projection of the particles distance
            // the origin is centered in the count position

            particles[partList[0]].active = false;
            particles[partList[0]].tau_exit = Tsim;
            particles[partList[0]].time  = Tsim;
            particles[partList[0]].tau_exitSampled = Tsim;

            particles[jPart].active = false;
            particles[jPart].tau_exit = Tsim;
            particles[jPart].time  = Tsim;
            particles[jPart].tau_exitSampled = Tsim;

            int tempList=0;
            for ( int n=j; n<N-1; n++ ){
                tempList=partList[n];
                partList[n]=partList[n+1];
                partList[n+1]=tempList;
            }
            for ( int n=0; n<N-1; n++ ){
                tempList=partList[n];
                partList[n]=partList[n+1];
                partList[n+1]=tempList;
            }


            break;
        }

    }

    particles[partList[0]].pos_exit[0] +=  gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt2TAU_BM;
    particles[partList[0]].pos_exit[1] +=  gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt2TAU_BM;
    particles[partList[0]].pos_exit[2] +=  gsl_ran_gaussian(r, 1) * particles[partList[0]].sqrtDiff * sqrt2TAU_BM;


    checkBound (particles[partList[0]].pos_exit, particles[partList[0]].pos_period, L );

    particles[partList[0]].tau_exit += tau_bm;

}

void BFstep_annih ( particle *particles, BFdistances *d, gsl_rng *r, double tau_bm, int N, double sqrt2TAU_BM, double L ) {
//dist,XYZ,deltaPos,varPos are just pointers to external free memory
//sqrtTAU_BM is sqrt(2*TAU_BM)
//deltaPos is an array of the increments in the position
//particles are in a box modeled with a soft core repulsion on the boundaries

    double dist,deltaPos [3], varPos[3];

    for (int i=0; i<N; i++) {

        particles[i].tau_exit += tau_bm;

        if (!particles[i].active)
            continue;

        deltaPos[0] = 0;
        deltaPos[1] = 0;
        deltaPos[2] = 0;


        for( int j=0; j<N; j++){

            if (!particles[j].active)
                continue;

            if (i==j) continue;
            // if the distance with another particle is lower than R_INTER we take into account their interaction
            if (i<j){
                dist = d[i].dd[j];
            }
            else{
                dist = d[j].dd[i];
            }

            if ( dist<particles[i].radius+particles[j].radius && particles[i].type==particles[j].type ){

                particles[i].active = false;
                particles[j].active = false;
                break;
            }

        }



        //brownian displacement
        deltaPos[0] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
        deltaPos[1] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;
        deltaPos[2] += gsl_ran_gaussian (r,1)*particles[i].sqrtDiff * sqrt2TAU_BM;

        particles[i].pos_exit[0] += deltaPos[0];
        particles[i].pos_exit[1] += deltaPos[1];
        particles[i].pos_exit[2] += deltaPos[2];
        checkBound (particles[i].pos_exit, particles[i].pos_period, L );


    }

}


