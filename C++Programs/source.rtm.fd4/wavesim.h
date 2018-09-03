#ifndef WAVESIM
#define WAVESIM
double **imcon(double **ug, double **dg, int nx, int ny, int nt);
double **ug_wve(double **rec,double *vel, int nx, int ny,double h, double dt,double T);
double **dg_wve(double *vel, int nx, int ny,int srcloc, double freq,double h, double dt,double T);
#endif
