#include <iostream>
#include <iomanip>
#include <accelmath.h>
#include "arrayman.h"
#include "hider.h"
#include "abcon.h"
#include "openacc.h"

#define idx(i,j,nx) (i*nx+j)

double **dg_wve(double *vel, int nx, int ny,int srcloc, double freq,double h, double dt,double T){
	int nt=int(T/dt);
	int Nx=nx+2;
	int Ny=ny+2;
	
	//Alloc array to store wavefield
	double **__restrict__ u=alloc_mat(nt,nx*ny);
	
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
	
	#pragma acc parallel loop copyin(stencil[0:nx*ny])
	for (int i=1;i<Ny-1;i++){
		for (int j=1;j<Nx-1;j++){
			stencil[(i-1)*nx+j-1]=i*Nx+j;
			}
		}

	//allocate array spatial operator for ABC
	int *__restrict__ left_sstep=new int[81];
	int *__restrict__ right_sstep=new int[81];
	int *__restrict__ bottom_sstep=new int[81];
	int *__restrict__ tstep=new int[81];

	//fill ABC operator
	gen_sstep(left_sstep,1);
	gen_sstep(right_sstep,-1);
	gen_sstep(bottom_sstep,-Nx);
	gen_tstep(tstep);

	//allocate array for coefficient of ABC hihgdon
	double **__restrict__ left_cfabc=alloc_mat(ny,81);
	double **__restrict__ right_cfabc=alloc_mat(ny,81);
	double **__restrict__ bottom_cfabc=alloc_mat(nx,81);
	
	#pragma acc kernels copyin(vel[0:nx*ny],left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81])
	for (int i=0;i<ny;i++)
		{
		gen_cfabc(left_cfabc[i],vel[i*nx],dt,h,beta);
		gen_cfabc(right_cfabc[i],vel[(i+1)*nx-1],dt,h,beta);
		}

	#pragma acc kernels copyin(vel[0:nx*ny],bottom_cfabc[0:nx][0:81])
	for (int i=0;i<nx;i++)
		{
		gen_cfabc(bottom_cfabc[i],vel[nx*ny-nx+i],dt,h,beta);
		}
	double t;
	
	#pragma acc enter data copyin(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],stencil[0:nx*ny],vel[0:nx*ny],left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],bottom_cfabc[0:nx][0:81],left_sstep[0:81],right_sstep[0:81],bottom_sstep[0:81],tstep[0:81])
	
	#pragma acc data copyout(u[0:nt][0:nx]) 
	for (int i=0; i<nt-1;i++)
	{
		//time step./w	
		t=dt*(i);
		if((i)%(nt/10)==0){
		std::cout<<std::fixed<<std::setprecision(1)<<"Calculating Wavefield ... "<<float(i)/float(nt)*100.0<<"%\n";}
		
		//Top neumann boundary
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny])
		for (int j=0;j<Nx;j++){
			U[3][j]=U[3][j+2*Nx];
			Ux[3][j]=Ux[3][j+2*Nx];
			Uy[3][j]=Uy[3][j+2*Nx];
				}
		
		#pragma acc kernels present(U[0:5][0:Nx*Ny])
		{
		int source_loc=stencil[srcloc];		
		U[3][source_loc]=-5.76*freq*freq*(1-16.0*(0.6*freq*t-1)*(0.6*freq*t-1)) *exp(-8.0* (0.6*freq*t-1)*(0.6*freq*t-1));
		}
		
		//Calculate Wavefield
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],stencil[0:nx*ny],vel[0:nx*ny])
		for (int j=0; j<nx*ny;j++){
			int pos=stencil[j];
			double cf1,cf2,cf3;
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			
			U[4][pos]=2.0*U[3][pos]-U[2][pos]+cf1*d2xd2y(U[3],h,pos,Nx)
						-cf2*(d4(U[3],Ux[3],h,pos,1)+d4(U[3],Uy[3],h,pos,Nx))
						+cf3*d2x2y(U[3],Ux[3],Uy[3],h,pos,Nx);
			
			Ux[4][pos]=2.0*Ux[3][pos]-Ux[2][pos]+cf1*d2xd2y(Ux[3],h,pos,Nx)
						-cf2*(d5(U[3],Ux[3],h,pos,1)+dx4y(U[3],Uy[3],h,pos,Nx))
						+cf3*d3x2y(U[3],Ux[3],h,pos,Nx);
			
			Uy[4][pos]=2.0*Uy[3][pos]-Uy[2][pos]+cf1*d2xd2y(Uy[3],h,pos,Nx)
						-cf2*(d4xy(U[3],Ux[3],h,pos,Nx)+d5(U[3],Uy[3],h,pos,Nx))
						+cf3*d2x3y(U[3],Uy[3],h,pos,Nx);
			}
		
		//calculate ABC higdon boundary
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],left_sstep[0:81],right_sstep[0:81],tstep[0:81])
		for (int j=1;j<Ny-1;j++) 
		{	
			U[4][j*Nx+1]=habc(U,left_cfabc[j-1],tstep,left_sstep,j*Nx+1);
			Ux[4][j*Nx+1]=habc(Ux,left_cfabc[j-1],tstep,left_sstep,j*Nx+1);
			Uy[4][j*Nx+1]=habc(Uy,left_cfabc[j-1],tstep,left_sstep,j*Nx+1);

			U[4][(j+1)*Nx-2]=habc(U,right_cfabc[j-1],tstep,right_sstep,(j+1)*Nx-2);
			Ux[4][(j+1)*Nx-2]=habc(Ux,right_cfabc[j-1],tstep,right_sstep,(j+1)*Nx-2);	
			Uy[4][(j+1)*Nx-2]=habc(Uy,right_cfabc[j-1],tstep,right_sstep,(j+1)*Nx-2);			

		}
		
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],tstep[0:81],bottom_cfabc[0:nx][0:81],bottom_sstep[0:81],tstep[0:81])
		for (int j=1;j<Nx-1;j++)
		{
			U[4][Ny*Nx-Nx+j]=habc(U,bottom_cfabc[j-1],tstep,bottom_sstep,Ny*Nx-Nx+j);
			Ux[4][Ny*Nx-Nx+j]=habc(U,bottom_cfabc[j-1],tstep,bottom_sstep,Ny*Nx-Nx+j);
			Uy[4][Ny*Nx-Nx+j]=habc(U,bottom_cfabc[j-1],tstep,bottom_sstep,Ny*Nx-Nx+j);
		}
		
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],stencil[0:nx*ny])
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
		for (int j=0;j<nx*ny;j++)
		{
			int pos=stencil[j];
			u[i+1][j]=U[4][pos];
		}
	}
	
	#pragma acc exit data delete(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny])
	#pragma acc exit data delete(left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],bottom_cfabc[0:nx][0:81])
	#pragma acc exit data delete(left_sstep[0:81],right_sstep[0:81],bottom_sstep[0:81],tstep[0:81])
	#pragma acc exit data delete(vel[0:nx*ny],stencil[0:nx*ny])
	
	free_mat_mem(U); free_mat_mem(Ux); free_mat_mem(Uy);
	free_mat_mem(left_cfabc); free_mat_mem(right_cfabc); free_mat_mem(bottom_cfabc);
	delete [] left_sstep; delete [] right_sstep; delete [] bottom_sstep; delete [] tstep;
	return u;}

double **ug_wve(double **rec,double *vel, int nx, int ny,double h,double dt,double T){
	int nt=int(T/dt);
	int Nx=nx+2;
	int Ny=ny+2;
	
	//Alloc array to store wavefield
	double **__restrict__ u=alloc_mat(nt,nx*ny);
	
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
	
	#pragma acc parallel loop copyin(stencil[0:nx*ny])
	for (int i=1;i<Ny-1;i++){
		for (int j=1;j<Nx-1;j++){
			stencil[(i-1)*nx+j-1]=i*Nx+j;
			}
		}

	//allocate array spatial operator for ABC
	int *__restrict__ left_sstep=new int[81];
	int *__restrict__ right_sstep=new int[81];
	int *__restrict__ bottom_sstep=new int[81];
	int *__restrict__ tstep=new int[81];

	//fill ABC operator
	gen_sstep(left_sstep,1);
	gen_sstep(right_sstep,-1);
	gen_sstep(bottom_sstep,-Nx);
	gen_tstep(tstep);

	//allocate array for coefficient of ABC hihgdon
	double **__restrict__ left_cfabc=alloc_mat(ny,81);
	double **__restrict__ right_cfabc=alloc_mat(ny,81);
	double **__restrict__ bottom_cfabc=alloc_mat(nx,81);
	
	#pragma acc kernels copyin(vel[0:nx*ny],left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81])
	for (int i=0;i<ny;i++)
		{
		gen_cfabc(left_cfabc[i],vel[i*nx],dt,h,beta);
		gen_cfabc(right_cfabc[i],vel[(i+1)*nx-1],dt,h,beta);
		}

	#pragma acc kernels copyin(vel[0:nx*ny],bottom_cfabc[0:nx][0:81])
	for (int i=0;i<nx;i++)
		{
		gen_cfabc(bottom_cfabc[i],vel[nx*ny-nx+i],dt,h,beta);
		}
	
	#pragma acc enter data copyin(rec[0:nt][0:nx*ny],U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],stencil[0:nx*ny],vel[0:nx*ny],left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],bottom_cfabc[0:nx][0:81],left_sstep[0:81],right_sstep[0:81],bottom_sstep[0:81],tstep[0:81])
	
	#pragma acc data copyout(u[0:nt][0:nx]) 
	for (int i=0; i<nt-1;i++)
	{
		//time step./w	
		t=dt*(i);
		if((i)%(nt/10)==0){
		std::cout<<std::fixed<<std::setprecision(1)<<"Calculating Wavefield ... "<<float(i)/float(nt)*100.0<<"%\n";}
		
		//Top neumann boundary
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny])
		for (int j=0;j<Nx;j++){
			U[3][j]=U[3][j+2*Nx];
			Ux[3][j]=Ux[3][j+2*Nx];
			Uy[3][j]=Uy[3][j+2*Nx];
				}
		
		#pragma acc kernels present(rec[0:nt][0:nx*ny],U[0:5][0:Nx*Ny])
		{
		for(int j=1;j<Nx-1;j++)
			{
				U[3][Nx+j]=rec[nt-1-i][j-1];
			}
		}
		
		//Calculate Wavefield
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],stencil[0:nx*ny],vel[0:nx*ny])
		for (int j=0; j<nx*ny;j++){
			int pos=stencil[j];
			double cf1,cf2,cf3;
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			
			U[4][pos]=2.0*U[3][pos]-U[2][pos]+cf1*d2xd2y(U[3],h,pos,Nx)
						-cf2*(d4(U[3],Ux[3],h,pos,1)+d4(U[3],Uy[3],h,pos,Nx))
						+cf3*d2x2y(U[3],Ux[3],Uy[3],h,pos,Nx);
			
			Ux[4][pos]=2.0*Ux[3][pos]-Ux[2][pos]+cf1*d2xd2y(Ux[3],h,pos,Nx)
						-cf2*(d5(U[3],Ux[3],h,pos,1)+dx4y(U[3],Uy[3],h,pos,Nx))
						+cf3*d3x2y(U[3],Ux[3],h,pos,Nx);
			
			Uy[4][pos]=2.0*Uy[3][pos]-Uy[2][pos]+cf1*d2xd2y(Uy[3],h,pos,Nx)
						-cf2*(d4xy(U[3],Ux[3],h,pos,Nx)+d5(U[3],Uy[3],h,pos,Nx))
						+cf3*d2x3y(U[3],Uy[3],h,pos,Nx);
			}
		
		//calculate ABC higdon boundary
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],left_sstep[0:81],right_sstep[0:81],tstep[0:81])
		for (int j=1;j<Ny-1;j++) 
		{	
			U[4][j*Nx+1]=habc(U,left_cfabc[j-1],tstep,left_sstep,j*Nx+1);
			Ux[4][j*Nx+1]=habc(Ux,left_cfabc[j-1],tstep,left_sstep,j*Nx+1);
			Uy[4][j*Nx+1]=habc(Uy,left_cfabc[j-1],tstep,left_sstep,j*Nx+1);

			U[4][(j+1)*Nx-2]=habc(U,right_cfabc[j-1],tstep,right_sstep,(j+1)*Nx-2);
			Ux[4][(j+1)*Nx-2]=habc(Ux,right_cfabc[j-1],tstep,right_sstep,(j+1)*Nx-2);	
			Uy[4][(j+1)*Nx-2]=habc(Uy,right_cfabc[j-1],tstep,right_sstep,(j+1)*Nx-2);			

		}
		
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],tstep[0:81],bottom_cfabc[0:nx][0:81],bottom_sstep[0:81],tstep[0:81])
		for (int j=1;j<Nx-1;j++)
		{
			U[4][Ny*Nx-Nx+j]=habc(U,bottom_cfabc[j-1],tstep,bottom_sstep,Ny*Nx-Nx+j);
			Ux[4][Ny*Nx-Nx+j]=habc(U,bottom_cfabc[j-1],tstep,bottom_sstep,Ny*Nx-Nx+j);
			Uy[4][Ny*Nx-Nx+j]=habc(U,bottom_cfabc[j-1],tstep,bottom_sstep,Ny*Nx-Nx+j);
		}
		
		#pragma acc parallel loop present(U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny],stencil[0:nx*ny])
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
		for (int j=0;j<nx*ny;j++)
		{
			int pos=stencil[j];
			u[i+1][j]=U[4][pos];
		}
	}
	
	#pragma acc exit data delete(rec[0:nt][0:nx*ny],U[0:5][0:Nx*Ny],Ux[0:5][0:Nx*Ny],Uy[0:5][0:Nx*Ny])
	#pragma acc exit data delete(left_cfabc[0:ny][0:81],right_cfabc[0:ny][0:81],bottom_cfabc[0:nx][0:81])
	#pragma acc exit data delete(left_sstep[0:81],right_sstep[0:81],bottom_sstep[0:81],tstep[0:81])
	#pragma acc exit data delete(vel[0:nx*ny],stencil[0:nx*ny])
	
	free_mat_mem(U); free_mat_mem(Ux); free_mat_mem(Uy);
	free_mat_mem(left_cfabc); free_mat_mem(right_cfabc); free_mat_mem(bottom_cfabc);
	delete [] left_sstep; delete [] right_sstep; delete [] bottom_sstep; delete [] tstep;
	return u;
}

double **imcon(double **ug, double **dg, int nx, int ny, int nt)
{
	double **__restrict__ I = alloc_mat(ny,nx);
	
	#pragma acc kernels copyin(I[0:ny][0:nx])
	for(int j=0;j<ny;j++)
		{
			for(int k=0;k<nx;k++)
			{
				I[j][k]=0.00;
			}
		}
	
	#pragma acc kernels copyin(ug[0:nt][0:nx*ny],dg[0:nt][0:nx*ny]) copyout(I[0:ny][0:nx])
	for (int i=0;i<nt;i++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int k=0;k<nx;k++)
			{
				int id=idx(j,k,nx);
				I[j][k]+=ug[nt-1-i][id]*dg[i][id];
			}
		}
	}
	return I;
}


