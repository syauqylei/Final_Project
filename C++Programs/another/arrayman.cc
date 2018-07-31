#include <bits/stdc++.h>
void pairsortd(int *a, double *b, int n)
{
    pair<int, double> pairt[n];
 
    // Storing the respective array
    // elements in pairs.
    for (int i = 0; i < n; i++) 
    {
        pairt[i].first = a[i];
        pairt[i].second = b[i];
    }
 
    // Sorting the pair array.
    sort(pairt, pairt + n);
     
    // Modifying original arrays
    for (int i = 0; i < n; i++) 
    {
        b[i] = pairt[i].second;
    }
}
void pairsort(int *a, int *b, int n)
{
    pair<int, int> pairt[n];
 
    // Storing the respective array
    // elements in pairs.
    for (int i = 0; i < n; i++) 
    {
        pairt[i].first = a[i];
        pairt[i].second = b[i];
    }
 
    // Sorting the pair array.
    sort(pairt, pairt + n);
     
    // Modifying original arrays
    for (int i = 0; i < n; i++) 
    {
        b[i] = pairt[i].second;
    }
}
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
