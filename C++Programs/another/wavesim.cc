#include <iostream>
#include <iomanip>
#include <cmath>
#include "arrayman.h"
#include "abcon.h"

double **wvenacd(double *vel, int nx, int ny,int srcloc, double freq,double h, double dt,double T){
	int nt=int(T/dt+1);
	int Nx=nx+2;
	int Ny=ny+2;
	
	std::cout<<"error?\n";
	double t;
	double **__restrict__ u=alloc_mat(nt,nx*ny);
	double **__restrict__ U=alloc_mat(nt+3,Nx*Ny);
	double **__restrict__ Ux=alloc_mat(nt+3,Nx*Ny);
	double **__restrict__ Uy=alloc_mat(nt+3,Nx*Ny);

	#pragma acc parallel loop 
	for (int i=0;i<3;i++){
		#pragma acc loop
		for(int j=0;j<Nx*Ny;j++){
			U[i][j]=0.0;
			Ux[i][j]=0.0;
			Uy[i][j]=0.0;
			}
		}
		//angle factor for abc higdon
	double *__restrict__ beta= new double[4];
	beta[0]=0.99;
	beta[1]=0.79;
	beta[2]=0.52;
	beta[3]=0.31;
	  std::cout<<"error?\n";

	int *__restrict__ stencil= new int[nx*ny];
	for (int i=1;i<Ny-1;i++){
		#pragma acc loop
		for (int j=1;j<Nx-1;j++){
			stencil[(i-1)*nx+j-1]=i*Nx+j;
			}
		}

	int *__restrict__ left_sstep=new int[81];
	int *__restrict__ right_sstep=new int[81];
	int *__restrict__ bottom_sstep=new int[81];
	int *__restrict__ tstep=new int[81];

	gen_sstep(left_sstep,1);
	gen_sstep(right_sstep,-1);
	gen_sstep(bottom_sstep,-Nx);
	gen_tstep(tstep);

	double **__restrict__ left_cfabc=alloc_mat(ny,81);
	double **__restrict__ right_cfabc=alloc_mat(ny,81);
	double **__restrict__ bottom_cfabc=alloc_mat(nx,81);
  std::cout<<"error?\n";
	
	for (int i=0;i<ny;i++)
		{	
		gen_cfabc(left_cfabc[i],vel[i*nx],dt,h,beta);
		gen_cfabc(right_cfabc[i],vel[(i+1)*nx-1],dt,h,beta);
		}
 
	for (int i=0;i<nx;i++)
		{
		gen_cfabc(bottom_cfabc[i],vel[nx*ny-nx+i],dt,h,beta);
		}
  std::cout<<"error?\n";

	//copyin (sstep,tstep,stencil
	for (int i=3; i<nt+2;i++){

  std::cout<<"error?\n";
		//time step./w	
		t=dt*(i-3);
		if((i-3)%(nt/10)==0){
		std::cout<<std::fixed<<std::setprecision(1)<<"Calculating Wavefield ... "<<float(i-3)/float(nt)*100.0<<"%\n";}
		//boundary Neumann 
		
		//independent
		#pragma acc parallel loop
		for (int j=0;j<Nx;j++){
			U[i][j]=U[i][j+2*Nx];
			Ux[i][j]=Ux[i][j+2*Nx];
			Uy[i][j]=Uy[i][j+2*Nx];
			}
		U[i][stencil[srcloc]]=-5.76*freq*freq*(1-16.0*(0.6*freq*t-1)*(0.6*freq*t-1)) *exp(-8.0* (0.6*freq*t-1)*(0.6*freq*t-1));
             
		//solution of nacd
		//independent wait upper loop 
		#pragma acc parallel loop
		for (int j=0; j<nx*ny;j++){
			int pos=stencil[j];
			double cf1,cf2,cf3;
			double D4x,D4y,D2x2y,UD2xD2y;
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			UD2xD2y=(U[i][pos+1]+U[i][pos-1]+U[i][pos+Nx]+U[i][pos-Nx]-4*U[i][pos])/h/h;
			D4x=-12.0/h/h/h/h*(U[i][pos+1]-2.0*U[i][pos]+U[i][pos-1])+(Ux[i][pos+1]-Ux[i][pos-1])*6.0/h/h/h;
			D4y=-12.0/h/h/h/h*(U[i][pos+Nx]-2.0*U[i][pos]+U[i][pos-Nx])+6.0/h/h/h*(Uy[i][pos+Nx]-Uy[i][pos-Nx]);
			D2x2y=1.0/h/h/h/h*(2.0*(U[i][pos+Nx]+U[i][pos+1]-2.0*U[i][pos]+U[i][pos-1]+U[i][pos-Nx])
					- U[i][pos-Nx+1]-U[i][pos-Nx-1]-U[i][pos+Nx+1]-U[i][pos+Nx-1])
					+1.0/2.0/h/h/h*(Ux[i][pos+Nx+1]+Ux[i][pos-Nx+1]-Ux[i][pos-Nx-1]-Ux[i][pos+Nx-1]-2.0*Ux[i][pos+1]+2.0*Ux[i][pos-1])
					+1.0/2.0/h/h/h*(Uy[i][pos+Nx+1]+Uy[i][pos+Nx-1]-Uy[i][pos-Nx-1]-Uy[i][pos-Nx+1]-2.0*Uy[i][pos+Nx]+2.0*Uy[i][pos-Nx]);
			
			U[i+1][pos]=2.0*U[i][pos]-U[i-1][pos]+cf1*UD2xD2y-cf2*(D4x+D4y)+cf3*D2x2y;
			u[i-3][j]=U[i+1][pos];//store wavefield
	
		}
		//independent
		#pragma acc parallel loop
		for(int j=0;j<nx*ny;j++){
			int pos=stencil[j];
			double cf1,cf2,cf3;
			double D5x,Dx4y,D3x2y,UxD2xD2y;
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			
			UxD2xD2y=(Ux[i][pos+1]+Ux[i][pos-1]+Ux[i][pos+Nx]+Ux[i][pos-Nx]-4*Ux[i][pos])/h/h;
			D5x=-90.0/h/h/h/h/h*(U[i][pos+1]-U[i][pos-1])+30.0/h/h/h/h*(Ux[i][pos+1]+4.0*Ux[i][pos]+Ux[i][pos-1]);
			Dx4y=-6.0/h/h/h/h/h*(U[i][pos+Nx+1]-U[i][pos-1-Nx]-U[i][pos-1+Nx]+U[i][pos+1-Nx]+2.0*U[i][pos-1]-2.0*U[i][pos+1])
					+3.0/h/h/h/h*(Uy[i][pos+Nx+1]+Uy[i][pos-Nx-1]-Uy[i][pos+Nx-1]-Uy[i][pos-Nx+1]);
			D3x2y=3.0/2.0/h/h/h/h/h*(U[i][pos+Nx+1]-U[i][pos-1-Nx]+U[i][pos+1-Nx]-U[i][pos-1+Nx]+2.0*U[i][pos-1]-2.0*U[i][pos+1])
					+3.0/2.0/h/h/h/h*(Ux[i][pos+Nx+1]+Ux[i][pos-Nx-1]+Ux[i][pos+Nx-1]+Ux[i][pos-Nx+1]-2.0*Ux[i][pos+1]-2.0*Ux[i][pos-1]);
			
			Ux[i+1][pos]=2.0*Ux[i][pos]-Ux[i-1][pos]+cf1*UxD2xD2y-cf2*(D5x+Dx4y)+cf3*D3x2y;
			}

		//independent
		#pragma acc parallel loop
		for(int j=0;j<nx*ny;j++){
			int pos=stencil[j];
			double D4xy,D5y,D2x3y,UyD2xD2y;
			double cf1,cf2,cf3;
			
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			
			UyD2xD2y=(Uy[i][pos+1]+Uy[i][pos-1]+Uy[i][pos+Nx]+Uy[i][pos-Nx]-4*Uy[i][pos])/h/h;
			D4xy=6.0/h/h/h/h/h*(U[i][pos+Nx+1]-U[i][pos-1-Nx]+U[i][pos-1+Nx]-U[i][pos+1-Nx]+2.0*U[i][pos-Nx]-2.0*U[i][pos+Nx])
					+3.0/h/h/h/h*(Ux[i][pos+Nx+1]+Ux[i][pos-Nx-1]-Ux[i][pos+Nx-1]-Ux[i][pos-Nx+1]);
			D5y=-90.0/h/h/h/h/h*(U[i][pos+Nx]-U[i][pos-Nx])+30.0/h/h/h/h*(Uy[i][pos+Nx]+4.0*Uy[i][pos]+Uy[i][pos-Nx]);
			D2x3y=-3.0/2.0/h/h/h/h/h*(U[i][pos+Nx+1]-U[i][pos-1-Nx]+U[i][pos-1+Nx]-U[i][pos+1-Nx]+2.0*U[i][pos-Nx]-2.0*U[i][pos+Nx])
					+3.0/2.0/h/h/h/h*(Uy[i][pos+Nx+1]+Uy[i][pos-Nx-1]+Uy[i][pos+Nx-1]+Uy[i][pos-Nx+1]-2.0*Uy[i][pos+Nx]-2.0*Uy[i][pos-Nx]);
					
			Uy[i+1][pos]=2.0*Uy[i][pos]-Uy[i-1][pos]+cf1*UyD2xD2y-cf2*(D4xy+D5y)+cf3*D2x3y;
			}
		//parallel}
		//independent wait
		#pragma acc parallel loop
		for (int j=1;j<Ny-1;j++)
		{	
			double Ubdrleft=0;
			#pragma acc loop
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=left_sstep[k];
				Ubdrleft+=-U[i+1+tshift][j*Nx+1+pos]*left_cfabc[j-1][k];}
				U[i+1][j*Nx+1]=Ubdrleft;

			double Uxbdrleft=0;
			#pragma acc loop 
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=left_sstep[k];
				Uxbdrleft+=-Ux[i+1+tshift][j*Nx+1+pos]*left_cfabc[j-1][k];}
				Ux[i+1][j*Nx+1]=Uxbdrleft;

			double Uybdrleft=0;
			#pragma acc loop
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=left_sstep[k];
				Uybdrleft+=-Uy[i+1+tshift][j*Nx+1+pos]*left_cfabc[j-1][k];}
				Uy[i+1][j*Nx+1]=Uybdrleft;

			double Ubdrright=0;
			#pragma acc loop
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Ubdrright+=-U[i+1+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];}
				U[i+1][(j-1)*Nx-2]=Ubdrright;			

			double Uxbdrright=0;
			#pragma acc loop
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Uxbdrright+=-Ux[i+1+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];}
				Ux[i+1][(j-1)*Nx-2]=Uxbdrright;			

			double Uybdrright=0;
			#pragma acc loop
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Uybdrright+=-Uy[i+1+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];}
				Uy[i+1][(j-1)*Nx-2]=Uybdrright;			

			}
		#pragma acc parallel loop
		for (int j=1;j<Nx-1;j++)
		{
			double Ubdrbottom=0;
			#pragma acc loop
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=bottom_sstep[k];
				Ubdrbottom+=-U[i+1+tshift][Ny*Nx-Nx+j+pos]*bottom_cfabc[j-1][k];}
				U[i+1][Ny*Nx-Nx+j]=Ubdrbottom;			
			
			double Uxbdrbottom=0;
			#pragma acc loop
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=bottom_sstep[k];
				Uxbdrbottom+=-Ux[i+1+tshift][Ny*Nx-Nx+j+pos]*bottom_cfabc[j-1][k];}
				Ux[i+1][Ny*Nx-Nx+j]=Uxbdrbottom;

			double Uybdrbottom=0;
			#pragma acc loop
			for (int k=0;k<81;k++)
				{
				int tshift=tstep[k];
				int pos=bottom_sstep[k];
				Uybdrbottom+=-Uy[i+1+tshift][Ny*Nx-Nx+j+pos]*bottom_cfabc[j-1][k];}
				Uy[i+1][Ny*Nx-Nx+j]=Uybdrbottom;
		
			}
		}
	free_mat_mem(left_cfabc);
	free_mat_mem(right_cfabc);
	free_mat_mem(bottom_cfabc);
	free_mat_mem(U);
	free_mat_mem(Ux);
	free_mat_mem(Uy);
	return u;
	}
