#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <ctime>
#include "rnw.h"
#include "arrayman.h"
#include "wavesim.h"

int main(int argc, char *argv[]){
	std::vector<double> vel;
	std::vector<int> srcloc;
	int nx,ny,nz,ns;
	double h,dt,T,freq;
	read_vel(argv[2],vel , nx,ny,nz,h);
	read_fwdset(argv[3],srcloc,ns,dt,T,freq);
	int nt=int(T/dt);
	double *Velocity=&vel[0];
	double **u;
	for (int i=0;i<ns;i++){
		std::cout<<"Create Record with Source Location ("<<srcloc[ns+i]<<","<<srcloc[i]<<")\n";
		int src_loc=srcloc[ns+i]*(nx)+srcloc[i];
		std::cout<<"Calculating Wave solution with FD2 .... \n";
		int start=std::clock();
		d=dg(Velocity,nx,ny,src_loc,freq,h,dt,T);
		u=
		int end=std::clock();
		std::cout << "time: " << (end-start)/double(CLOCKS_PER_SEC)*1000 << std::endl;
		write_rec_txt(argv[1], u,srcloc[i],srcloc[ns+i],h,dt, nx,ny,nt);
		}
	free_mat_mem(u);
	}
