#include <iostream>
#include <vector>
#include "rnw.h"
#include "arrayman.h"
int main (int argc,char *argv[])
{
	std::vector<double> rec;
	int nx,nt,xsource,ysource;
	double h,dt;
	std::string ext="pr";
	read_rec(argv[1],rec);
	std::string fname=argv[1]+ext;
	read_rec_pr(fname,nx,nt,xsource,ysource,h,dt);
	double **record=alloc_mat(nt,nx);
	double v_water=1400.0;
	rec_2d(record,rec,nx,nt);

	int ntl=(xsource*h/v_water)/dt+(11*h/v_water)/dt;
	if (ntl>nt){ntl=nt;}
	for (int i=0;i<ntl;i++)
	{
		double t=i*dt;
		int x_step=v_water*t/h;
		int xl=xsource+1+10-x_step;
		for (int j=0;j<xl;j++)
		{
			record[i][j]=0.00;
		}
	}
	int ntr=((nx-xsource)*h/v_water)/dt+(11*h/v_water)/dt;
	if (ntr>nt){ntr=nt;}
	for (int i=0;i<ntr;i++)
	{
		double t=i*dt;
		int x_step=v_water*t/h;
		int xr=xsource-10+x_step;
		
		for (int j=xr;j<nx;j++)
		{
			record[i][j]=0.00;
		}
	}
	
	std::string filename=argv[2];
	write_rec_txt(filename,record,xsource,ysource,h,dt,nx,nt);
	free_mat_mem(record);
	return 0;
}
