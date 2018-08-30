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
	std::vector<double> rec;
	
	int nx,ny,nz,ns;
	double h,dt,T,freq;
	
	read_vel(argv[2],vel , nx,ny,nz,h);
	read_fwdset(argv[3],srcloc,ns,dt,T,freq);
	read_rec(argv[4],rec);
	
	int nt=int(T/dt);
	double **record=alloc_mat(nt,nx);
	
	rec_2d(record,rec,nx,nt);
	
	int nt=int(T/dt);
	double *Velocity=&vel[0];
	double **u;
	double **d;
	double **I;
	for (int i=0;i<ns;i++){
		std::cout<<"Create Image from Record with Source Location ("<<srcloc[ns+i]<<","<<srcloc[i]<<")\n";
		
		int src_loc=srcloc[ns+i]*(nx)+srcloc[i];
		
		std::cout<<"Calculating Wave solution with NACD .... \n";
		
		int start=std::clock();
		
		d=dg_wve(Velocity,nx,ny,src_loc,freq,h,dt,T);
		u=ug_wve(record,Velocity,nx,ny,h,dt,T);
		I=imcon(u,d,nx,ny,nt);
		
		int end=std::clock();
		
		std::cout << "time: " << (end-start)/double(CLOCKS_PER_SEC)*1000 << std::endl;
		write_img_txt(argv[1],I,nx,ny);
		}
	free_mat_mem(u);
	free_mat_mem(d);
	free_mat_mem(I);
	
	}
