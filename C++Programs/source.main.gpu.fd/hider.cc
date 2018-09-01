#pragma acc routine seq
double d2xd2y(double *U, double h, int pos, int Nx)
{
	return (U[pos+1]+U[pos-1]+U[pos+Nx]+U[pos-Nx]-4*U[pos])/h/h;
}

#pragma acc routine seq
double d4(double *U, double *U_, double h, int pos,int dom)
{
	return -12.0/h/h/h/h*(U[pos+dom]-2.0*U[pos]+U[pos-dom])+(U_[pos+dom]-U_[pos-dom])*6.0/h/h/h;
}

#pragma acc routine seq
double d5(double *U, double *U_, double h, int pos,int dom)
{
	return -90.0/h/h/h/h/h*(U[pos+dom]-U[pos-dom])+30.0/h/h/h/h*(U_[pos+dom]+4.0*U_[pos]+U_[pos-dom]);
}

#pragma acc routine seq
double d2x2y(double *U, double *Ux, double *Uy, double h, int pos, int Nx)
{
return 1.0/h/h/h/h*(2.0*(U[pos+Nx]+U[pos+1]-2.0*U[pos]+U[pos-1]+U[pos-Nx])
		- U[pos-Nx+1]-U[pos-Nx-1]-U[pos+Nx+1]-U[pos+Nx-1])
		+1.0/2.0/h/h/h*(Ux[pos+Nx+1]+Ux[pos-Nx+1]-Ux[pos-Nx-1]-Ux[pos+Nx-1]-2.0*Ux[pos+1]+2.0*Ux[pos-1])
		+1.0/2.0/h/h/h*(Uy[pos+Nx+1]+Uy[pos+Nx-1]-Uy[pos-Nx-1]-Uy[pos-Nx+1]-2.0*Uy[pos+Nx]+2.0*Uy[pos-Nx]);
}

#pragma acc routine seq
double dx4y(double *U, double *Uy, double h, int pos, int Nx)
{
return -6.0/h/h/h/h/h*(U[pos+Nx+1]-U[pos-1-Nx]-U[pos-1+Nx]+U[pos+1-Nx]+2.0*U[pos-1]-2.0*U[pos+1])
		+3.0/h/h/h/h*(Uy[pos+Nx+1]+Uy[pos-Nx-1]-Uy[pos+Nx-1]-Uy[pos-Nx+1]);
}

#pragma acc routine seq
double d4xy(double *U, double *Ux, double h, int pos, int Nx)
{
return -6.0/h/h/h/h/h*(U[pos+Nx+1]-U[pos-1-Nx]+U[pos-1+Nx]-U[pos+1-Nx]+2.0*U[pos-Nx]-2.0*U[pos+Nx])
		+3.0/h/h/h/h*(Ux[pos+Nx+1]+Ux[pos-Nx-1]-Ux[pos+Nx-1]-Ux[pos-Nx+1]);
}

#pragma acc routine seq
double d3x2y(double *U, double *Ux, double h, int pos, int Nx)
{
	return 3.0/2.0/h/h/h/h/h*(U[pos+Nx+1]-U[pos-1-Nx]+U[pos+1-Nx]-U[pos-1+Nx]+2.0*U[pos-1]-2.0*U[pos+1])
					+3.0/2.0/h/h/h/h*(Ux[pos+Nx+1]+Ux[pos-Nx-1]+Ux[pos+Nx-1]+Ux[pos-Nx+1]-2.0*Ux[pos+1]-2.0*Ux[pos-1]);
}

#pragma acc routine seq
double d2x3y(double *U, double *Uy, double h, int pos, int Nx)
{
	return-3.0/2.0/h/h/h/h/h*(U[pos+Nx+1]-U[pos-1-Nx]+U[pos-1+Nx]-U[pos+1-Nx]+2.0*U[pos-Nx]-2.0*U[pos+Nx])
					+3.0/2.0/h/h/h/h*(Uy[pos+Nx+1]+Uy[pos-Nx-1]+Uy[pos+Nx-1]+Uy[pos-Nx+1]-2.0*Uy[pos+Nx]-2.0*Uy[pos-Nx]);
}
