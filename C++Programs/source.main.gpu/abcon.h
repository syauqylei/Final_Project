#ifndef ABC
#define ABC
#pragma acc routine
void gen_sstep(int *s_step,int pole);
#pragma acc routine
void gen_tstep(int *t_step);
#pragma acc routine
void gen_cfabc(double *cfabc,double c,double dt,double h,double *beta);
#endif
