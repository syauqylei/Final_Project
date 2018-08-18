#include <bits/stdc++.h>
double **alloc_mat(const int nrows, const int ncols){
	double **mat=new double*[nrows];
	
	mat[0]= new double[nrows*ncols];
	for (int i=1; i<nrows;i++){
		mat[i]=&mat[0][i*ncols];
	}	
	return mat;
}

void free_mat_mem(double **mat){
	delete [] mat[0];
	delete [] mat;
}
