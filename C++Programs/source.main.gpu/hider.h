#ifndef HIDER
#define	HIDER

#pragma acc routine seq
extern double d2xd2y(double *U, double h, int pos, int Nx);

#pragma acc routine seq
extern double d4(double *U, double *U_, double h, int pos, int dom);

#pragma acc routine seq
extern double d5(double *U, double *U_, double h, int pos, int dom);

#pragma acc routine seq
extern double d2x2y(double *U, double *Ux, double *Uy, double h, int pos, int Nx);

#pragma acc routine seq
extern double dx4y(double *U, double *Uy, double h, int pos, int Nx);

#pragma acc routine seq
extern double d4xy(double *U, double *Ux, double h, int pos, int Nx);

#pragma acc routine seq
extern double d3x2y(double *U, double *Ux, double h, int pos, int Nx);

#pragma acc routine seq
extern double d2x3y(double *U, double *Uy, double h, int pos, int Nx)

#endif
