// author luigisbailo


#pragma once

void run_annih_PQ ( int N_A, int N_B, int R_A, int R_B, double D_A, double D_B, double tau_bm, double alpha, double Tsim, int nProj, double L,  int *stat, int * n_annih ) {

	const gsl_rng_type *Type;
	gsl_rng *r;
	gsl_rng_env_setup ();
	Type = gsl_rng_default;
	r = gsl_rng_alloc (Type);
	FILE *devurandom = fopen("/dev/urandom","r");
	unsigned long int seed;
	fread(&seed, sizeof(seed), 1, devurandom);
	fclose(devurandom);
	gsl_rng_set(r, seed);

	double tProj = Tsim / nProj;
	int countProj = 1;

	int N = N_A + N_B;

	double R;

	const double sqrt2TAU_BM = sqrt(2*tau_bm);

	double Teq = tProj;

	double distRow [N], maxSh;

	particle particles [N];
	double shells [N];
	int partList [N];

	initPos_hybGF ( particles, r, N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, L );

	initShell_GF ( particles, r, N, tau_bm, sqrt2TAU_BM, L, &stat[1]);

	//sort() is a prebuild c++ funct. It sorts particles for increasing exit times
	std::sort ( particles, particles+N, compareTime );
	for (int n=0; n<N; n++) partList[n]=n;


	while ( particles[partList[0]].tau_exit < Tsim ) {


		if ( particles[partList[0]].burst )
			stat[0]++;

		updatePart_GF ( &particles[partList[0]], r, tau_bm, L );
		//differently from aGF, updatePart() here samples also the exit position from the shell

		getDist ( particles, partList, distRow, &maxSh, N, L );

		burst_PQ_GF ( particles, partList, distRow, r, N, partList[0], L);

		if (particles[partList[0]].tau_exitSampled<particles[partList[0]].time) {
			R = getR_GF(particles, partList, shells, distRow, N, L);
		}
		else{
			R=0;
		}

		particles[partList[0]].burst = false; //when a particle is burst its position is updated and put on top of the list

		if ( R > 0 ) {

			stat [1] ++;
			if (R>L/20) R=L/20;
			GFstep_GF ( &particles[partList[0]], r, R );
			particles[partList[0]].gf = true;

		}
		else{

			stat [2] ++;
			BMstep_annih ( particles, partList, distRow, r, tau_bm,  sqrt2TAU_BM, N, L, Tsim );
			particles[partList[0]].gf = false;
		}

		sortBurst ( particles, partList, N);

		//The particles are sorted according tau_exit, while instead tau_exitSampled keeps only memory of a previously sampled exit time
		sortPart ( particles,partList,N);


		if ( particles[partList[0]].tau_exit > tProj | particles[partList[0]].tau_exit == tProj ) {


			synchPart_PQ_GF ( particles, partList, r, N, tProj, L );


			int n_active = 0;

			for ( int n=0; n<N; n++){
				if (particles[n].active)
					n_active++;
			}

			n_annih [countProj-1] = n_active;

			initShell_GF ( particles, r, N, tau_bm, sqrt2TAU_BM, L, &stat[1]);

			std::sort ( particles, particles+N, compareTime );
			for (int n=0; n<N; n++) partList[n]=n;

			countProj ++;
			tProj = countProj * Tsim / nProj;

		}

	} ;

	gsl_rng_free (r);

}
