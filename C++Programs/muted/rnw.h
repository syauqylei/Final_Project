#ifndef RNW
#define RNW
void read_rec(const std::string& filename, std::vector <double>& rec);
void read_rec_pr(const std::string& filename, int &nx,int &nt,int &xsource,int &ysource, double &h,double &dt);
void rec_2d(double **rec_out,std::vector <double>& rec, int nx,int nt);
void write_rec_txt(const std::string& filename, double **U,int xsource,int ysource,double h,double dt , int nx, int nt);
#endif
