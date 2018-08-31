#ifndef RNW
#define RNW
void read_vel(const std::string& filename, std::vector <double>& vel , int &nx,int &ny, int &nz, double &h);
void w_dat(const std::string& filename, double *Vel_mod, double **P,double dt,double h,int nt,int nx, int ny, int nz,double a,double b, double c);
void read_fwdset(const std::string& filename, std::vector <int>& srcloc,int &ns, double &dt,double &T,double &fm);
void write_wve_txt(const std::string& filename, double **U , int nx,int ny, int nt);
void write_rec_txt(const std::string& filename, double **U,int xsource,int ysource,double h,double dt , int nx,int ny, int nt);
#endif
