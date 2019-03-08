// author luigisbailo


void fig1 () {

	double b=1;
	double D=1;
	double tau = 1;
 	double t1 = 0.01;
 	double t2 = 0.1;
 	double t3 = 0.99;
	double Pder_t1, Pder_t2, Pder_t3, PQ00der_t1, PQ00der_t2, PQ00der_t3;

	for ( double r=b/100; r<b; r+=b/100){

		Pder_t1 = Pder(r,t1,b,D,Sfunct(t1,b,D));
		Pder_t2 = Pder(r,t2,b,D,Sfunct(t2,b,D));
		Pder_t3 = Pder(r,t3,b,D,Sfunct(t3,b,D));
		PQ00der_t1 = PQ00der (r,t1,tau,b,D,Sder(tau,b,D));
		PQ00der_t2 = PQ00der (r,t2,tau,b,D,Sder(tau,b,D));
		PQ00der_t3 = PQ00der (r,t3,tau,b,D,Sder(tau,b,D));

		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n", r, Pder_t1, Pder_t2, Pder_t3, PQ00der_t1, PQ00der_t2, PQ00der_t3 );

	}

}
