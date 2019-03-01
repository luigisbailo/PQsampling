// author luigisbailo


#pragma once

void printPos ( particle *particles, int *partList, int N){

  std::cout <<"Part.\t"<<"Time\t"<<"Tau\t"<<"Shell\t"
       <<"x\t"<<"y\t"<<"z\t"
       <<"x_ex\t"<<"y_ex\t"<<"z_ex\t"
       <<"Radius\t" <<"BURST\n";

  for (int count=0; count<N; count++){
    std::cout <<particles[partList[count]].label << "\t"
	 <<particles[partList[count]].time << "\t" 
	 <<particles[partList[count]].tau_exit<<"\t"
	 <<particles[partList[count]].shell<<"\t"
	 <<particles[partList[count]].pos[0]<<"\t"
	 <<particles[partList[count]].pos[1]<<"\t"
	 <<particles[partList[count]].pos[2] <<"\t"
	 <<particles[partList[count]].pos_exit[0]<<"\t"
	 <<particles[partList[count]].pos_exit[1]<<"\t"
	 <<particles[partList[count]].pos_exit[2] <<"\t"
	 <<particles[partList[count]].radius <<"\t"
	 <<particles[partList[count]].burst <<"\n";
  }
}


void printPosInit ( particle *particles, int *partList, int N){

    std::cout <<"Part.\t"<<"Time\t"<<"Tau\t"<<"Shell\t"
              <<"x\t"<<"y\t"<<"z\t"
              <<"x_ex\t"<<"y_ex\t"<<"z_ex\t"
              <<"Radius\t" <<"BURST\n";

    for (int count=0; count<N; count++){
        std::cout <<particles[partList[count]].label << "\t"
                  <<particles[partList[count]].time << "\t"
                  <<particles[partList[count]].tau_exit<<"\t"
                  <<particles[partList[count]].shell<<"\t"
                  <<particles[partList[count]].pos[0]<<"\t"
                  <<particles[partList[count]].pos[1]<<"\t"
                  <<particles[partList[count]].pos[2] <<"\t"
                  <<particles[partList[count]].pos_exit[0]<<"\t"
                  <<particles[partList[count]].pos_exit[1]<<"\t"
                  <<particles[partList[count]].pos_exit[2] <<"\t"
                  <<particles[partList[count]].radius <<"\t"
                  <<particles[partList[count]].burst <<"\n";
    }
}


void printPosPQ ( particle *particles, int *partList, int N){

	std::cout <<"Part.\t"<<"Time\t"<<"Tau\t"<<"Shell\t" <<"BURST\n";

	for (int count=0; count<N; count++){
		std::cout <<particles[partList[count]].label << "\t"
				  <<particles[partList[count]].time << "\t"
				  <<particles[partList[count]].tau_exit<<"\t"
 				  <<particles[partList[count]].tau_exitSampled<<"\t"
				  <<particles[partList[count]].shell<<"\t"
				  <<particles[partList[count]].radius <<"\t"
				  <<particles[partList[count]].burst <<"\n";
	}
}

void printDist_per (particle *particles, int *partList, int N, double L ) {
 
  for (int n=1; n<N; n++){
	double dist = sqrt(dist2_per(&particles[partList[0]],&particles[partList[n]],L));
	double distNext = sqrt(dist2next_per(&particles[partList[0]],&particles[partList[n]],L));
	double ishell = particles[partList[0]].shell;
	double jshell = particles[partList[n]].shell;
	std::cout << "R "<<particles[partList[0]].label<< "-" <<particles[partList[n]].label<<" = "<< dist << "\t"
	     << "R - shells -inter = "<<dist - (ishell+jshell+particles[partList[n]].radius+particles[partList[0]].radius) << "\t"
	     << "Rnext -inter = "<<distNext - (particles[partList[n]].radius+particles[partList[0]].radius) 
	     << "\n";  

  }

}



void print_shell ( particle *particles, int N){
  
  for (int i=0; i<N; i++){
    std::cout << "Particle " << i << "\tDiff = " << particles[i].Diff
	 << "\tR_bm = " << particles[i].R_bd
	 << "\tShell = " << particles[i].shell 
	 << "\tTau = " << particles[i].tau_exit << "\n";
            
  }
}


