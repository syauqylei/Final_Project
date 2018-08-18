#include <iostream>
#include <iomanip>
#include <cmath>

const double b = 0.0;
double qx(double bta,double cfv){
	double rslt=(b*(bta+cfv)-cfv)/((bta+cfv)*(1-b));
	return rslt;
	}

double qt(double bta,double cfv){
	double rslt=(b*(bta+cfv)-bta)/((bta+cfv)*(1-b));
	return rslt;
	}

void gen_sstep(int *s_step,int pole){
	int itr[3]={3,9,27};
	s_step[0]=0;
	s_step[1]=pole*1;
	s_step[2]=0;
	
	for (int i=0;i<3;i++){
		for (int j=0; j<2;j++){
			for (int k=0; k<itr[i];k++){
				s_step[itr[i]+j*itr[i]+k]=(s_step[k]+s_step[j+1]);
				}
			}
		}
	}

void gen_tstep(int *t_step){
	int itr[3]={3,9,27};
	t_step[0]=0;
	t_step[1]=0;
	t_step[2]=-1;

	for (int i=0;i<3;i++){
		for (int j=0; j<2;j++){
			for (int k=0; k<itr[i];k++){
				t_step[itr[i]+j*itr[i]+k]=t_step[k]+t_step[j+1];
				}
			}
		}
	}

#pragma acc routine seq
void gen_cfabc(double *cfabc,double c,double dt,double h,double *beta){
	double cfl=c*dt/h;
	double abc_eq[12];
	int iter[3]={3,9,27};
	for (int i=0;i<4;i++){
			abc_eq[i*3]=1;
			abc_eq[i*3+1]=qx(beta[i],cfl);
			abc_eq[i*3+2]=qt(beta[i],cfl);
			}
	cfabc[0]=abc_eq[0];
	cfabc[1]=abc_eq[1];
	cfabc[2]=abc_eq[2];
	for (int i=0;i<3;i++){
		int iterate=iter[i];
		for (int j=0; j<2;j++){
			int iterliml=iter[i];
			int iterlimr=j*iter[i];
			for (int k=0; k<iterate;k++){
				cfabc[iterliml+iterlimr+k]=cfabc[k]*abc_eq[i*3+3+j+1];
				}
			}
		}
	cfabc[0]=0;
	}

#pragma acc routine seq
double habc(double **U, double *cfabc, int *tstep, int *sstep, int pos)
{
	double Ubdr=0;
	#pragma acc loop reduction(+:Ubdr)
	for (int k=0;k<81;k++)
	{
		int tshift=tstep[k];
		int sshift=sstep[k];
		Ubdr+=-U[4+tshift][pos+sshift]*cfabc[k];
	}
	return Ubdr;
}
