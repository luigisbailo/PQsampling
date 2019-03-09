// author luigisbailo


void fig2_annih (){

    double D_A = 0.01;
    double D_B = 0.01;
    double R_A = 2.5;
    double R_B = 2.5;

    double tau_bm = 0.1;
    const int N = 10;
    const int N_A = 5;
    const int N_B = 5;

    double alpha= 9;
    double L = 14;

    int nProj=20;
    double Tsim=1000;
    int Nsamples=100;

    int stat[3];
    int n_stat[nProj];
    double n_stat_P [Nsamples][nProj];
    double n_stat_PQ [Nsamples][nProj];
    double n_stat_BM [Nsamples][nProj];


    for ( int count = 0; count < Nsamples; count++ ){


        run_annih_P ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, nProj, L, stat, n_stat );

        for ( int t=0; t<nProj; t++){

            n_stat_P [count][t] = n_stat[t];

        }


        run_annih_PQ ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, alpha, Tsim, nProj, L, stat, n_stat );

        for ( int t=0; t<nProj; t++){

            n_stat_PQ [count][t] = n_stat[t];

        }


        run_annih_BM ( N_A, N_B, R_A, R_B, D_A, D_B, tau_bm, Tsim, nProj, L, n_stat );

        for ( int t=0; t<nProj; t++){

            n_stat_BM [count][t] = n_stat[t];
        }


    }


    double av_n_P [nProj];
    double av_n_PQ [nProj];
    double av_n_BM [nProj];

    for ( int t=0; t<nProj; t++ ){

        av_n_P [t] = 0;
        av_n_PQ [t] = 0;
        av_n_BM [t] = 0;

        for ( int count=0; count<Nsamples; count++){

            av_n_P [t] += n_stat_P [count][t];
            av_n_PQ [t] += n_stat_PQ [count][t];
            av_n_BM [t] += n_stat_BM [count][t];

        }

        av_n_P [t] = av_n_P [t] / Nsamples;
        av_n_PQ [t] = av_n_PQ [t] / Nsamples;
        av_n_BM [t] = av_n_BM [t] / Nsamples;


    }


    double sd_n_P [nProj];
    double sd_n_PQ [nProj];
    double sd_n_BM [nProj];

    for ( int t=0; t<nProj; t++ ){

        sd_n_P [t] = 0;
        sd_n_PQ [t] = 0;
        sd_n_BM [t] = 0;

        for ( int count=0; count<Nsamples; count++){

            sd_n_P [t] += pow(n_stat_P[count][t]-av_n_P[t],2);
            sd_n_PQ [t] += pow(n_stat_PQ[count][t]-av_n_PQ[t],2);
            sd_n_BM [t] += pow(n_stat_BM[count][t]-av_n_BM[t],2);

        }

        sd_n_P [t] = sqrt(sd_n_P [t] / pow(Nsamples,2) );
        sd_n_PQ [t] = sqrt(sd_n_PQ [t] / pow(Nsamples,2) );
        sd_n_BM [t] = sqrt(sd_n_BM [t] / pow(Nsamples,2) );

    }



    for ( int t=0; t<nProj; t++){

        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n",
               (t+1)*Tsim/nProj,  av_n_P[t], av_n_PQ[t], av_n_BM[t], sd_n_P[t], sd_n_PQ[t], sd_n_BM[t]);


    }

}

