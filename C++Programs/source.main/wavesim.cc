#include <iostream>
#include <iomanip>
#include <cmath>
#include "arrayman.h"
#include "abcon.h"

#define index(i,j,N) (i*N+j)

double **wvenacd(double *vel, int nx, int ny,int srcloc, double freq,double h, double dt,double T){
	int nt=int(T/dt);
	int Nx=nx+2;
	int Ny=ny+2;
	//Alloc array to store wavefield
	double **__restrict__ u=alloc_mat(nt,nx);

	double **__restrict__ U=alloc_mat(5,Nx*Ny);
	double **__restrict__ Ux=alloc_mat(5,Nx*Ny);
	double **__restrict__ Uy=alloc_mat(5,Nx*Ny);

	//allocate array for Un-1 Un Un-1
	double *__restrict__ Uo= new double[Nx*Ny];
	double *__restrict__ Up= new double[Nx*Ny];
	double *__restrict__ Um= new double[Nx*Ny];

	//allocate array for Uxn-1 Uxn Uxn-1
	double *__restrict__ Uxo= new double[Nx*Ny];
	double *__restrict__ Uxp= new double[Nx*Ny];
	double *__restrict__ Uxm= new double[Nx*Ny];
	
	//allocate array for Uyn-1 Uyn Uyn-1
	double *__restrict__ Uyo= new double[Nx*Ny];
	double *__restrict__ Uyp= new double[Nx*Ny];
	double *__restrict__ Uym= new double[Nx*Ny];
	//initial condition
	
	#pragma acc kernels copyin(Uxp[:Ny*Nx],Uym[:Ny*Nx],Uyo[:Ny*Nx],Uyp[:Ny*Nx],Up[:Ny*Nx],Uxm[:Ny*Nx],Uxo[:Ny*Nx],Um[:Ny*Nx],Uo[:Ny*Nx])
	for(int j=0;j<Nx*Ny;j++)
	{
		Up[j]=0.0;
		Uo[j]=0.0;
		Um[j]=0.0;
		Uxp[j]=0.0;
		Uxo[j]=0.0;
		Uxm[j]=0.0;
		Uyp[j]=0.0;
		Uyo[j]=0.0;
		Uym[j]=0.0;
	}
	
	//create angle factor for abc higdon
	double *__restrict__ beta= new double[4];
	beta[0]=0.99;
	beta[1]=0.79;
	beta[2]=0.52;
	beta[3]=0.31;
	
	//create stencil
	int *__restrict__ stencil= new int[nx*ny];
	#pragma acc kernels copyin(stencil[:nx*ny])
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

	//fill coefficient to array
	#pragma acc kernels
	for (int i=0;i<ny;i++)
		{	
		gen_cfabc(left_cfabc[i],vel[i*nx],dt,h,beta);
		gen_cfabc(right_cfabc[i],vel[(i+1)*nx-1],dt,h,beta);
		}
	#pragma acc kernels
	for (int i=0;i<nx;i++)
		{
		gen_cfabc(bottom_cfabc[i],vel[nx*ny-nx+i],dt,h,beta);
		}
	
	double t;
	
	#pragma acc data copyin(left_sstep[:81],right_sstep[:81],bottom_sstep[:81],tstep[:81],left_cfabc[:ny][:81],right_cfabc[:ny][:81],bottom_cfabc[:nx][:81])		
	for (int i=0; i<nt-1;i++)
	{
		//time step./w	
		t=dt*(i);
		if((i)%(nt/10)==0){
		std::cout<<std::fixed<<std::setprecision(1)<<"Calculating Wavefield ... "<<float(i)/float(nt)*100.0<<"%\n";}
		
		
		//Top neuman boundary
		#pragma acc parallel loop
		for (int j=0;j<Nx;j++){
			Uo[j]=Uo[j+2*Nx];
			Uxo[j]=Uxo[j+2*Nx];
			Uyo[j]=Uyo[j+2*Nx];
			}
		
		//Source
		int source_loc=stencil[srcloc];
		Uo[source_loc]=-5.76*freq*freq*(1-16.0*(0.6*freq*t-1)*(0.6*freq*t-1)) *exp(-8.0* (0.6*freq*t-1)*(0.6*freq*t-1));
        
        //Calculate Wavefield
		#pragma acc parallel loop 
		for (int j=0; j<nx*ny;j++){
			int pos=stencil[j];
			double cf1,cf2,cf3;
			double D4x,D4y,D2x2y,UD2xD2y;
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			UD2xD2y=(Uo[pos+1]+Uo[pos-1]+Uo[pos+Nx]+Uo[pos-Nx]-4*Uo[pos])/h/h;
			D4x=-12.0/h/h/h/h*(Uo[pos+1]-2.0*Uo[pos]+Uo[pos-1])+(Uxo[pos+1]-Uxo[pos-1])*6.0/h/h/h;
			D4y=-12.0/h/h/h/h*(Uo[pos+Nx]-2.0*Uo[pos]+Uo[pos-Nx])+6.0/h/h/h*(Uyo[pos+Nx]-Uyo[pos-Nx]);
			D2x2y=1.0/h/h/h/h*(2.0*(Uo[pos+Nx]+Uo[pos+1]-2.0*Uo[pos]+Uo[pos-1]+Uo[pos-Nx])
					- Uo[pos-Nx+1]-Uo[pos-Nx-1]-Uo[pos+Nx+1]-Uo[pos+Nx-1])
					+1.0/2.0/h/h/h*(Uxo[pos+Nx+1]+Uxo[pos-Nx+1]-Uxo[pos-Nx-1]-Uxo[pos+Nx-1]-2.0*Uxo[pos+1]+2.0*Uxo[pos-1])
					+1.0/2.0/h/h/h*(Uyo[pos+Nx+1]+Uyo[pos+Nx-1]-Uyo[pos-Nx-1]-Uyo[pos-Nx+1]-2.0*Uyo[pos+Nx]+2.0*Uyo[pos-Nx]);
			Up[pos]=2.0*Uo[pos]-Um[pos]+cf1*UD2xD2y-cf2*(D4x+D4y)+cf3*D2x2y;
		}
		
		#pragma acc parallel loop 
		for(int j=0;j<nx*ny;j++){
			int pos=stencil[j];
			double cf1,cf2,cf3;
			double D5x,Dx4y,D3x2y,UxD2xD2y;
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			UxD2xD2y=(Uxo[pos+1]+Uxo[pos-1]+Uxo[pos+Nx]+Uxo[pos-Nx]-4*Uxo[pos])/h/h;
			D5x=-90.0/h/h/h/h/h*(Uo[pos+1]-Uo[pos-1])+30.0/h/h/h/h*(Uxo[pos+1]+4.0*Uxo[pos]+Uxo[pos-1]);
			Dx4y=-6.0/h/h/h/h/h*(Uo[pos+Nx+1]-Uo[pos-1-Nx]-Uo[pos-1+Nx]+Uo[pos+1-Nx]+2.0*Uo[pos-1]-2.0*Uo[pos+1])
					+3.0/h/h/h/h*(Uyo[pos+Nx+1]+Uyo[pos-Nx-1]-Uyo[pos+Nx-1]-Uyo[pos-Nx+1]);
			D3x2y=3.0/2.0/h/h/h/h/h*(Uo[pos+Nx+1]-Uo[pos-1-Nx]+Uo[pos+1-Nx]-Uo[pos-1+Nx]+2.0*Uo[pos-1]-2.0*Uo[pos+1])
					+3.0/2.0/h/h/h/h*(Uxo[pos+Nx+1]+Uxo[pos-Nx-1]+Uxo[pos+Nx-1]+Uxo[pos-Nx+1]-2.0*Uxo[pos+1]-2.0*Uxo[pos-1]);
			Uxp[pos]=2.0*Uxo[pos]-Uxm[pos]+cf1*UxD2xD2y-cf2*(D5x+Dx4y)+cf3*D3x2y;
			}
		
		#pragma acc parallel loop
		for(int j=0;j<nx*ny;j++){
			int pos=stencil[j];
			double D4xy,D5y,D2x3y,UyD2xD2y;
			double cf1,cf2,cf3;
			
			cf1=vel[j]*vel[j]*dt*dt;
			cf2=(vel[j]*vel[j]*dt*dt*h*h-vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/12.0;
			cf3=(vel[j]*vel[j]*vel[j]*vel[j]*dt*dt*dt*dt)/6.0;
			
			UyD2xD2y=(Uyo[pos+1]+Uyo[pos-1]+Uyo[pos+Nx]+Uyo[pos-Nx]-4*Uyo[pos])/h/h;
			D4xy=6.0/h/h/h/h/h*(Uo[pos+Nx+1]-Uo[pos-1-Nx]+Uo[pos-1+Nx]-Uo[pos+1-Nx]+2.0*Uo[pos-Nx]-2.0*Uo[pos+Nx])
					+3.0/h/h/h/h*(Uxo[pos+Nx+1]+Uxo[pos-Nx-1]-Uxo[pos+Nx-1]-Uxo[pos-Nx+1]);
			D5y=-90.0/h/h/h/h/h*(Uo[pos+Nx]-Uo[pos-Nx])+30.0/h/h/h/h*(Uyo[pos+Nx]+4.0*Uyo[pos]+Uyo[pos-Nx]);
			D2x3y=-3.0/2.0/h/h/h/h/h*(Uo[pos+Nx+1]-Uo[pos-1-Nx]+Uo[pos-1+Nx]-Uo[pos+1-Nx]+2.0*Uo[pos-Nx]-2.0*Uo[pos+Nx])
					+3.0/2.0/h/h/h/h*(Uyo[pos+Nx+1]+Uyo[pos-Nx-1]+Uyo[pos+Nx-1]+Uyo[pos-Nx+1]-2.0*Uyo[pos+Nx]-2.0*Uyo[pos-Nx]);
					
			Uyp[pos]=2.0*Uyo[pos]-Uym[pos]+cf1*UyD2xD2y-cf2*(D4xy+D5y)+cf3*D2x3y;
			}
		
		//store wavefield for ABC calculation
		#pragma acc data copyin(U[:5][:Nx*ny],Ux[:5][:Nx*ny],Uy[:5][:Nx*ny])
		{
		#pragma acc parallel loop present(U[:5][:Nx*ny],Ux[:5][:Nx*ny],Uy[:5][:Nx*ny])
		for(int j=0; j<Nx*Ny;j++)
		{
			U[0][j]=U[1][j];
			U[1][j]=U[2][j];
			U[2][j]=Um[j];
			U[3][j]=Uo[j];
			U[4][j]=Up[j];
			
			Ux[0][j]=Ux[1][j];
			Ux[1][j]=Ux[2][j];
			Ux[2][j]=Uxm[j];
			Ux[3][j]=Uxo[j];
			Ux[4][j]=Uxp[j];
		
			Uy[0][j]=Uy[1][j];
			Uy[1][j]=Uy[2][j];
			Uy[2][j]=Uym[j];
			Uy[3][j]=Uyo[j];
			Uy[4][j]=Uyp[j];
		}

		//calculate ABC higdon boundary
		
		#pragma acc parallel loop present(left_sstep[:81],right_sstep[:81],tstep[:81],left_cfabc[:ny][:81],right_cfabc[:ny][:81])
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
				Up[j*Nx+1]=Ubdrleft;

			double Uxbdrleft=0;
			#pragma acc loop reduction(+:Uxbdrleft)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=left_sstep[k];
				Uxbdrleft+=-Ux[4+tshift][j*Nx+1+pos]*left_cfabc[j-1][k];
			}
				Uxp[j*Nx+1]=Uxbdrleft;

			double Uybdrleft=0;
			#pragma acc loop reduction(+:Uybdrleft)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=left_sstep[k];
				Uybdrleft+=-Uy[4+tshift][j*Nx+1+pos]*left_cfabc[j-1][k];
			}
				Uyp[j*Nx+1]=Uybdrleft;

			double Ubdrright=0;
			#pragma acc loop reduction(+:Ubdrright)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Ubdrright+=-U[4+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];
			}
				Up[(j-1)*Nx-2]=Ubdrright;			

			double Uxbdrright=0;
			#pragma acc loop reduction(+:Uxbdrright)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Uxbdrright+=-Ux[4+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];
			}
				Uxp[(j-1)*Nx-2]=Uxbdrright;			

			double Uybdrright=0;
			#pragma acc loop reduction(+:Uybdrright)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=right_sstep[k];
				Uybdrright+=-Uy[4+tshift][(j-1)*Nx-2+pos]*right_cfabc[j-1][k];
			}
				Uyp[(j-1)*Nx-2]=Uybdrright;			

		}
		
		#pragma acc parallel loop present(bottom_sstep[:81],tstep[:81],bottom_cfabc[:nx][:81]) 
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
				Up[Ny*Nx-Nx+j]=Ubdrbottom;			
			
			double Uxbdrbottom=0;
			#pragma acc loop reduction(+:Uxbdrbottom)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=bottom_sstep[k];
				Uxbdrbottom+=-Ux[4+tshift][Ny*Nx-Nx+j+pos]*bottom_cfabc[j-1][k];
			}
				Uxp[Ny*Nx-Nx+j]=Uxbdrbottom;

			double Uybdrbottom=0;
			#pragma acc loop reduction(+:Uybdrbottom)
			for (int k=0;k<81;k++)
			{
				int tshift=tstep[k];
				int pos=bottom_sstep[k];
				Uybdrbottom+=-Uy[4+tshift][Ny*Nx-Nx+j+pos]*bottom_cfabc[j-1][k];
			}
				Uyp[Ny*Nx-Nx+j]=Uybdrbottom;
		}
		}
		
		#pragma acc parallel loop 
		for (int j=0;j<nx*ny;j++)	
		{
			int pos=stencil[j];
			Um[pos]=Uo[pos];
			Uo[pos]=Up[pos];
			Uxm[pos]=Uxo[pos];
			Uxo[pos]=Uxp[pos];
			Uym[pos]=Uyo[pos];
			Uyo[pos]=Uyp[pos];
		}
		#pragma acc kernels
		for (int j=0;j<nx;j++)
		{
			int pos=stencil[j];
			u[i+1][j]=Up[pos];
		}
	}
	#pragma acc exit data delete(U,Ux,Uy,left_sstep,right_sstep,bottom_sstep,tstep,left_cfabc,right_cfabc,bottom_cfabc,Up,Uo,Um,Uxo,Uxm,Uxp,Uyo,Uyp,Uym)
	free_mat_mem(U); free_mat_mem(Ux); free_mat_mem(Uy);
	free_mat_mem(left_cfabc); free_mat_mem(right_cfabc); free_mat_mem(bottom_cfabc);
	delete [] left_sstep; delete [] right_sstep; delete [] bottom_sstep; delete [] tstep;
	delete [] Up; delete [] Uo; delete [] Um;
	delete [] Uxp; delete [] Uxo; delete [] Uxm;
	delete [] Uyp; delete [] Uyo; delete [] Uym;

	return u;}
