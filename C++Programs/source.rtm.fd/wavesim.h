#ifndef WAVESIM
#define WAVESIM
double **wvenacd(double *vel, int nx, int ny,int srcloc, double fmax,double h, double dt,double T);
double **wvefd(double *vel, int nx, int ny,int srcloc, double fmax,double h, double dt,double T);
#endif
