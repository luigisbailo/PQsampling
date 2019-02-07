#pragma once

double getR_GF ( particle *particles, int *particleList, double *shells, double *distRow, int N, double L ){

  int iPart,jPart;
  double Rshell;

  iPart = particleList[0];

  //shells[0] contains the minimial shell
  if ( particles[iPart].gf == true )  
    shells[0] = particles[iPart].R_bd; 
  else 
    shells[0] = particles[iPart].R_gfrd;

  for ( int j=1; j<N; j++) {

    jPart = particleList [j];
    // Rshell is an array of 3 elements, the 1st-[0] contains the distance with the j-particle, the 2nd the distance with the expected exit position, the 3rd
    if ( particles[jPart].gf == false ){
      Rshell =  distRow [j] / 2;
    }
    else {
      Rshell = distRow [j] - particles[jPart].shell;
    }

    if( Rshell < shells[0] )
      return 0;
    else 
      shells[j] = Rshell;

  }

  return *std::min_element( shells+1,shells+N);

}


