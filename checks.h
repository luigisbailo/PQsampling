// author luigisbailo


#pragma once

void check_aGF (  struct particle *particles, int *partList, int N, double L ){
  

  if (particles[partList[0]].time > particles[partList[0]].tau_exit){
    printf("ERROR: exit time\n");
    exit (EXIT_FAILURE);
  }

  if ( particles[partList[0]].pos[0]>L | particles[partList[0]].pos[0]<0 | particles[partList[0]].pos[1]>L | particles[partList[0]].pos[1] <0  | particles[partList[0]].pos[2]>L | particles[partList[0]].pos[2]<0  ) {
    printf("ERROR: pos\n");
    exit (EXIT_FAILURE);
  
  } 

  if (  particles[partList[0]].gf == 0  &&
       (particles[partList[0]].pos_exit[0]>L | particles[partList[0]].pos_exit[0]<0 | particles[partList[0]].pos_exit[1]>L | particles[partList[0]].pos_exit[1] <0  | particles[partList[0]].pos_exit[2]>L | particles[partList[0]].pos_exit[2]<0 ) ) {
    printf("ERROR: pos_exit\n");
    exit (EXIT_FAILURE);
  
  } 

  for (int count=1; count<N; count++){
    int jPart = partList[count];
    if ( sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell - particles[partList[0]].radius - particles[jPart].radius < - 0.0000000001  && particles[jPart].shell>0 ){
      printf("ERROR: distance particles");
      exit(EXIT_FAILURE);
    }
  }

}


void check_GF (  struct particle *particles, int *partList, int N, double L ){
  

  if (particles[partList[0]].time > particles[partList[0]].tau_exit){
    printf("ERROR: exit time\n");
    exit (EXIT_FAILURE);
  }


  if ( particles[partList[0]].pos[0]>L | particles[partList[0]].pos[0]<0 | particles[partList[0]].pos[1]>L | particles[partList[0]].pos[1] <0  | particles[partList[0]].pos[2]>L | particles[partList[0]].pos[2]<0  ) {
    printf("ERROR: pos\n");
    exit (EXIT_FAILURE);
  
  } 

  for (int count=1; count<N; count++){

    int jPart = partList[count];
    if ( sqrt(dist2_per(&particles[partList[0]],&particles[jPart],L)) - particles[partList[0]].shell - particles[jPart].shell - particles[partList[0]].radius - particles[jPart].radius < - 0.0000000001  && particles[jPart].shell>0 ){
      printf("ERROR: distance particles\n");
      exit(EXIT_FAILURE);
    }

  }
  
}


void check_times ( struct particle *particles, int *partList, int N) {

  for (int count=1; count<N; count++){

    if ( particles[count].tau_exit<particles[partList[0]].time){
      printf("ERROR: time1\n");
      exit (EXIT_FAILURE);
    }
    if ( particles[partList[0]].time<particles[count].time){
      printf("ERROR: time2\n");
      exit (EXIT_FAILURE);
    }

  }

}


