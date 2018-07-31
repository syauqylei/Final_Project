#ifndef ARRAYMAN_H
#define ARRAYMAN_H
double **alloc_mat(const int nrows, const int ncols);
void free_mat_mem(double **mat);
void pairsortd(int *a,double *b, int n);
void pairsort(int *a,int *b, int n);
#endif

