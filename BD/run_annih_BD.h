#pragma once

void run_annih_BM ( int N_A, int N_B, int R_A, int R_B, double D_A, double D_B, double tau_bm, double Tsim, int nProj, double L, int *n_annih) {


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

	const double sqrt2TAU_BM = sqrt(2*tau_bm);

	const int N = N_A + N_B;

	struct BFdistances d[N];
	struct particle particles [N];

	double tProj = Tsim / nProj;
	int countProj = 1; 

    initPos_BM ( particles, r, N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, L); 

	do {

		getBFdistances ( particles, d, N ,L );

    	BFstep_annih ( particles, d, r, tau_bm, N, sqrt2TAU_BM, L );

		BFupdate ( particles, N);

		if ( particles[0].time > tProj ) {

			int n_active = 0;

			for ( int n=0; n<N; n++){
				if (particles[n].active == 0)
					n_active++;
			}

			n_annih [countProj-1] = n_active;

			countProj ++;
			tProj = countProj * Tsim / nProj;

		}

	} while ( particles[0].time < Tsim );


    gsl_rng_free (r);


}
