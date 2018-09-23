#pragma acc routine seq
double d2xd2y(double *U, double h, int pos, int Nx, int Ny)
{
	int top_idc;
	int bot_idc;
	
	if ( pos-2*Nx < 0 ) 
	{
		top_idc=pos+2*Nx;
	}
	
	else 
	{
		top_idc=pos-2*Nx;
	}
	
	if ( pos+2*Nx > Nx*Ny-Nx-1 ) 
	{
		bot_idc=pos+Nx;
	}
	
	else 
	{
		bot_idc=pos+2*Nx;
	}
	
	return (U[pos+2]+U[pos-2]+U[bot_idc]+U[top_idc]+16.0*U[pos+1]+16.0*U[pos-1]+16.0*U[pos+Nx]+16.0*U[pos-Nx]-60.0*U[pos])/h/h/12.0;
}
