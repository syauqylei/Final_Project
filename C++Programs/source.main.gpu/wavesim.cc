#include <iostream>
#include <iomanip>
#include <accelmath.h>
#include "arrayman.h"
#include "hider.h"
#include "abcon.h"
#include "openacc.h"

double **wvenacd(double *vel, int nx, int ny,int srcloc, double freq,double h, double dt,double T){
	int nt=int(T/dt);
	int Nx=nx+2;
	int Ny=ny+2;
	
	#pragma acc enter data copyin(vel[0:nx*ny])
	//Alloc array to store wavefield
	double **__restrict__ u=alloc_mat(nt,nx);
	double **__restrict__ U=alloc_mat(5,Nx*Ny);
	double **__restrict__ Ux=alloc_mat(5,Nx*Ny);
	double **__restrict__ Uy=alloc_mat(5,Nx*Ny);
	
	//create data on device
	#pragma acc kernels copyout(U[:5][:Nx*Ny],Ux[:5][:Nx*Ny],Uy[:5][:Nx*Ny])
	for (int i=0;i<5;i++)
	{
		for (int j=0;j<Nx*Ny;j++)
		{
			U[i][j]=0.0;
			Ux[i][j]=0.0;
			Uy[i][j]=0.0;
		}
	}
	
	
	//create angle factor for abc higdon
	double *__restrict__ beta= new double[4];
	beta[0]=0.99;
	beta[1]=0.79;
	beta[2]=0.52;
	beta[3]=0.31;
	
	//create stencil
	int *__restrict__ stencil= new int[nx*ny];
	#pragma acc enter data create(stencil[0:nx*ny])
	
	for (int i=1;i<Ny-1;i++){
		for (int j=1;j<Nx-1;j++){
			stencil[(i-1)*nx+j-1]=i*Nx+j;
			}
		}

	#pragma acc update device(stencil[0:nx*ny])
	//allocate array spatial operator for ABC
	int *__restrict__ left_sstep=new int[81];
	int *__restrict__ right_sstep=new int[81];
	int *__restrict__ bottom_sstep=new int[81];
	int *__restrict__ tstep=new int[81];

	#pragma acc enter data create(left_sstep[0:81],right_sstep[0:81],bottom_sstep[0:81],tstep[0:81])
	//fill ABC operator
	gen_sstep(left_sstep,1);
	gen_sstep(right_sstep,-1);
	gen_sstep(bottom_sstep,-Nx);
	gen_tstep(tstep);

	#pragma acc update device(left_sstep[0:81],right_sstep[0:81],bottom_sstep[0:81],tstep[0:81])
	
	//allocate array for coefficient of ABC hihgdon
	double **__restrict__ left_cfabc=alloc_mat(ny,81);
	double **__restrict__ right_cfabc=alloc_mat(ny,81);
	double **__restrict__ bottom_cfabc=alloc_mat(nx,81);
	
	#pragma acc enter data create(left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],bottom_cfabc[0:nx][0:81])
	
	#pragma acc kernels copyin(vel[0:nx*ny]) present(left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81])
	for (int i=0;i<ny;i++)
		{
		gen_cfabc(left_cfabc[i],vel[i*nx],dt,h,beta);
		gen_cfabc(right_cfabc[i],vel[(i+1)*nx-1],dt,h,beta);
		}

	#pragma acc kernels copyin(vel[0:nx*ny]) present(bottom_cfabc[0:nx][0:81])
	for (int i=0;i<nx;i++)
		{
		gen_cfabc(bottom_cfabc[i],vel[nx*ny-nx+i],dt,h,beta);
		}
	double t;
	#pragma acc data copyin(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny]) copyout(u[0:nt][0:nx]) present(stencil[0:nx*ny],vel[0:nx*ny],left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],bottom_cfabc[0:nx][0:81],left_sstep[0:81],right_sstep[0:81],bottom_sstep[0:81],tstep[0:81])
	for (int i=0; i<nt-1;i++)
	{
		//time step./w	
		t=dt*(i);
		if((i)%(nt/10)==0){
		std::cout<<std::fixed<<std::setprecision(1)<<"Calculating Wavefield ... "<<float(i)/float(nt)*100.0<<"%\n";}
		
		//Top neumann boundary
		#pragma acc parallel loop
		for (int j=0;j<Nx;j++){
			U[3][j]=U[3][j+2*Nx];
			Ux[3][j]=Ux[3][j+2*Nx];
			Uy[3][j]=Uy[3][j+2*Nx];
				}
		
		#pragma acc kernels
		{
		int source_loc=stencil[srcloc];		
		U[3][source_loc]=-5.76*freq*freq*(1-16.0*(0.6*freq*t-1)*(0.6*freq*t-1)) *exp(-8.0* (0.6*freq*t-1)*(0.6*freq*t-1));
		}
		
		//Calculate Wavefield
		#pragma acc parallel loop 
		for (int j=0; j<nx*ny;j++){
			int pos=stencil[j];
			double cf1,cf2,cf3;
			double D4x,D4y,D2x2y,UD2xD2y;
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			
			UD2xD2y=d2xd2y(U[3],h,pos,Nx);
			D4x=d4(U[3],Ux[3],h,pos,1);
			D4y=d4(U[3],Uy[3],h,pos,Nx);
			D2x2y=d2x2y(U[3],Ux[3],Uy[3],h,pos,Nx);
			U[4][pos]=2.0*U[3][pos]-U[2][pos]+cf1*UD2xD2y-cf2*(D4x+D4y)+cf3*D2x2y;
		}
		
		#pragma acc parallel loop
		for(int j=0;j<nx*ny;j++){
			int pos=stencil[j];
			double cf1,cf2,cf3;
			double D5x,Dx4y,D3x2y,UxD2xD2y;
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			
			UxD2xD2y=d2xd2y(Ux[3],h,pos,Nx);
			D5x=d5(U[3],Ux[3],h,pos,1);
			Dx4y=dx4y(U[3],Uy[3],h,pos,Nx);
			D3x2y=d3x2y(U[3],Ux[3],h,pos,Nx);
			Ux[4][pos]=2.0*Ux[3][pos]-Ux[2][pos]+cf1*UxD2xD2y-cf2*(D5x+Dx4y)+cf3*D3x2y;
			}
		
		#pragma acc parallel loop
		for(int j=0;j<nx*ny;j++){
			int pos=stencil[j];
			double D4xy,D5y,D2x3y,UyD2xD2y;
			double cf1,cf2,cf3;
			
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			
			UyD2xD2y=d2xd2y(Uy[3],h,pos,Nx);
			D4xy=d4xy(U[3],Ux[3],h,pos,Nx);
			D5y=d5(U[3],Uy[3],h,pos,Nx);
			D2x3y=d2x3y(U[3],Uy[3],h,pos,Nx);
					
			Uy[4][pos]=2.0*Uy[3][pos]-Uy[2][pos]+cf1*UyD2xD2y-cf2*(D4xy+D5y)+cf3*D2x3y;
			}
		
		//calculate ABC higdon boundary
		#pragma acc parallel loop
		for (int j=1;j<Ny-1;j++)
		{	
			double Ubdrleft=0;
			#pragma acc loop reduction(+:Ubdrleft)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=left_sstep[k];
				Ubdrleft+=-U[4+tshift][j*Nx+1+pos]*left_cfabc[j-1][k];
			}
				U[4][j*Nx+1]=Ubdrleft;

			double Uxbdrleft=0;
			#pragma acc loop reduction(+:Uxbdrleft)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=left_sstep[k];
				Uxbdrleft+=-Ux[4+tshift][j*Nx+1+pos]*left_cfabc[j-1][k];
			}
				Ux[4][j*Nx+1]=Uxbdrleft;

			double Uybdrleft=0;
			#pragma acc loop reduction(+:Uybdrleft)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=left_sstep[k];
				Uybdrleft+=-Uy[4+tshift][j*Nx+1+pos]*left_cfabc[j-1][k];
			}
				Uy[4][j*Nx+1]=Uybdrleft;

			double Ubdrright=0;
			#pragma acc loop reduction(+:Ubdrright)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Ubdrright+=-U[4+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];
			}
				U[4][(j+1)*Nx-2]=Ubdrright;			

			double Uxbdrright=0;
			#pragma acc loop reduction(+:Uxbdrright)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Uxbdrright+=-Ux[4+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];
			}
				Ux[4][(j+1)*Nx-2]=Uxbdrright;			

			double Uybdrright=0;
			#pragma acc loop reduction(+:Uybdrright)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Uybdrright+=-Uy[4+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];
			}
				Uy[4][(j+1)*Nx-2]=Uybdrright;			

		}
		
		#pragma acc parallel loop
		for (int j=1;j<Nx-1;j++)
		{
			double Ubdrbottom=0;
			#pragma acc loop reduction(+:Ubdrbottom)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=bottom_sstep[k];
				Ubdrbottom+=-U[4+tshift][Ny*Nx-Nx+j+pos]*bottom_cfabc[j-1][k];
			}
				U[4][Ny*Nx-Nx+j]=Ubdrbottom;			
			
			double Uxbdrbottom=0;
			#pragma acc loop reduction(+:Uxbdrbottom)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=bottom_sstep[k];
				Uxbdrbottom+=-Ux[4+tshift][Ny*Nx-Nx+j+pos]*bottom_cfabc[j-1][k];
			}
				Ux[4][Ny*Nx-Nx+j]=Uxbdrbottom;

			double Uybdrbottom=0;
			#pragma acc loop reduction(+:Uybdrbottom)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=bottom_sstep[k];
				Uybdrbottom+=-Uy[4+tshift][Ny*Nx-Nx+j+pos]*bottom_cfabc[j-1][k];
			}
				Uy[4][Ny*Nx-Nx+j]=Uybdrbottom;
		}
		
		#pragma acc parallel loop
		for (int j=0;j<nx*ny;j++)	
		{
			int pos=stencil[j];
			U[0][pos]=U[1][pos];
			U[1][pos]=U[2][pos];
			U[2][pos]=U[3][pos];
			U[3][pos]=U[4][pos];
			
			Ux[0][pos]=Ux[1][pos];
			Ux[1][pos]=Ux[2][pos];
			Ux[2][pos]=Ux[3][pos];
			Ux[3][pos]=Ux[4][pos];
			
			Uy[0][pos]=Uy[1][pos];
			Uy[1][pos]=Uy[2][pos];
			Uy[2][pos]=Uy[3][pos];
			Uy[3][pos]=Uy[4][pos];
		}
		
		#pragma acc parallel loop
		for (int j=0;j<nx;j++)
		{
			int pos=stencil[j];
			u[i+1][j]=U[4][pos];
		}
	}
	
	#pragma acc exit data delete(left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],bottom_cfabc[0:nx][0:81])
	#pragma acc exit data delete(left_sstep[0:81],right_sstep[0:81],bottom_sstep[0:81],tstep[0:81])
	#pragma acc exit data delete(vel[0:nx*ny])
	free_mat_mem(U); free_mat_mem(Ux); free_mat_mem(Uy);
	free_mat_mem(left_cfabc); free_mat_mem(right_cfabc); free_mat_mem(bottom_cfabc);
	delete [] left_sstep; delete [] right_sstep; delete [] bottom_sstep; delete [] tstep;
	return u;}
