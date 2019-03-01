// author luigisbailo


#pragma once

void check_aGF (  particle *particles, int *partList, int N, double L ){
  

  if (particles[partList[0]].time > particles[partList[0]].tau_exit){
    std::cout << "ERROR: exit time\n";
    exit (EXIT_FAILURE);
  }

  if ( particles[partList[0]].pos[0]>L | particles[partList[0]].pos[0]<0 | particles[partList[0]].pos[1]>L | particles[partList[0]].pos[1] <0  | particles[partList[0]].pos[2]>L | particles[partList[0]].pos[2]<0  ) {
    std::cout << "ERROR: pos\n";
    exit (EXIT_FAILURE);
  
  } 

  if (  particles[partList[0]].gf == true  &&
       (particles[partList[0]].pos_exit[0]>L | particles[partList[0]].pos_exit[0]<0 | particles[partList[0]].pos_exit[1]>L | particles[partList[0]].pos_exit[1] <0  | particles[partList[0]].pos_exit[2]>L | particles[partList[0]].pos_exit[2]<0 ) ) {
    std::cout << "ERROR: pos_exit\n";
    exit (EXIT_FAILURE);
  
  } 

  for (int count=1; count<N; count++){
    int jPart = partList[count];
    if ( sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell - particles[partList[0]].radius - particles[jPart].radius < - 0.0000000001  && particles[jPart].shell>0 ){
      std::cout << "ERROR: distance particles " << particles[partList[0]].label << "-" <<  particles[jPart].label << "\n";
      std::cout << sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell << "\n";
      exit(EXIT_FAILURE);
    }
  }

}


void check_GF (  particle *particles, int *partList, int N, double L ){
  

  if (particles[partList[0]].time > particles[partList[0]].tau_exit){
    std::cout << "ERROR: exit time\n";
    exit (EXIT_FAILURE);
  }


  if ( particles[partList[0]].pos[0]>L | particles[partList[0]].pos[0]<0 | particles[partList[0]].pos[1]>L | particles[partList[0]].pos[1] <0  | particles[partList[0]].pos[2]>L | particles[partList[0]].pos[2]<0  ) {
    std::cout << "ERROR: pos\n";
    // printPos (particles,N);
    exit (EXIT_FAILURE);
  
  } 

  for (int count=1; count<N; count++){

    int jPart = partList[count];
    if ( sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell - particles[partList[0]].radius - particles[jPart].radius < - 0.0000000001  && particles[jPart].shell>0 ){
      std::cout << "ERROR: distance particles " << particles[partList[0]].label << "-" <<  particles[jPart].label << "\n";
      std::cout << sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell << "\n";
      exit(EXIT_FAILURE);
    }

  }
  
}


void check_times ( particle *particles, int *partList, int N) {

  for (int count=1; count<N; count++){

    if ( particles[count].tau_exit<particles[partList[0]].time){
      std::cout << "ERROR: time1\n";
      std::cout << particles[partList[0]].time << "\t" << particles[count].tau_exit << "\n";
      exit (EXIT_FAILURE);
    }
    if ( particles[partList[0]].time<particles[count].time){
      std::cout << "ERROR: time2\n";
      std::cout << particles[partList[0]].time << "\t" << particles[count].time << "\n";
    printPos (particles,partList,N);
      exit (EXIT_FAILURE);
    }

  }

}


