#ifndef ABC
#define ABC
void gen_sstep(int *s_step,int pole);
void gen_tstep(int *t_step);

#pragma acc routine seq
extern void gen_cfabc(double *cfabc,double c,double dt,double h,double *beta);

#endif
