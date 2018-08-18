#ifndef ABC
#define ABC

#pragma acc routine seq
extern void gen_sstep(int *s_step,int pole);

#pragma acc routine seq
extern void gen_tstep(int *t_step);

#pragma acc routine seq
extern void gen_cfabc(double *cfabc,double c,double dt,double h,double *beta);
#endif
