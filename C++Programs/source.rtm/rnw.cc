#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "arrayman.h"

void read_vel(const std::string& filename, std::vector <double>& vel , int &nx,int &ny, int &nz, double &h){
	std::ifstream f(filename);
	
	int i=0;
	double content;
	while (f >> content)
	{
		i++;
		if(i==1){
			nx=content;}
		else if(i==2){
			ny=content;}
		else if(i==3){
			nz=content;}
		else if(i==4){
			h=content;}
		else{vel.push_back(content);}
	}
	}


void w_dat(const std::string& filename, double *Vel_mod, double **P,double dt,double h,int nt,int nx, int ny, int nz,double a,double b, double c){
	std::string ext=".dat";
	std::string fname;
	std::ofstream myfile;
	int tprint=(nt/10);
	std::cout << "Start to write File .....\n";
	for (int t=0; t<nt;t++){
		if (t%tprint==0){
		std::cout << "Writing File ..... "<<(float)t/nt*100 <<"%\n";}
		std::stringstream ss;
		ss <<std::setfill('0')<< std::setw(5) <<t;
		fname=filename+ss.str()+ ext;
		myfile.open(fname);		
		myfile << "TITLE = "<<"\""<<fname<<"\""<<std::endl;
		myfile << "VARIABLES = "<<"\"x\""<<","<<"\"y\""<<","<<"\"z\""<<","<<"\"Wavefield\""<<","<<"\"Velocity\""<<std::endl;
		myfile << "ZONE T = "<<"\"Frame "<<std::to_string(t)<<"\""<<", I = "<< std::to_string(nx)<<", J = "<< std::to_string(ny)<<", K = "<< std::to_string(nz)<<std::endl;
		for (int i=0; i<nz;i++){
			for (int j=0; j<ny;j++){
				for (int k=0;k<nx;k++){
					myfile << std::fixed<< std::setprecision(5)<<a*k*h+h/2.0<<","<<b*j*h+h/2.0<<","<<c*i*h+h/2.0<<","<<P[t][i*nx*nz+j*nx+k]<<","<<Vel_mod[i*nx*nz+j*nx+k]<<"\n";
				}
			}
		}
	myfile.close();
	}
}


void read_fwdset(const std::string& filename, std::vector <int>& srcloc,int &ns, double &dt,double &T,double &fm){
	std::ifstream f(filename);
	std::vector<double> source_loc;
	int i=0;
	double content;while (f >> content)
	{
		i++;
		if(i==1){
			dt=content;}
		else if(i==2){
			T=content;}
		else if(i==3){
			fm=content;}
		else{srcloc.push_back(content);}
	}
	ns=srcloc.size()/2;
	}

void write_wve_txt(const std::string& filename, double **U , int nx,int ny, int nt){
	std::string ext=".txt";
	std::string fname=filename+ext;
	std::ofstream file;
	file.open(fname);
	
	for(int i=0;i<nt;i++){
		if(i%(nt/10)==0){
		std::cout<<"Writing Wavefield ... "<<float(i)/float(nt)*100<<"%\n";}
		for(int j=0;j<ny;j++){
			for(int k=0;k<nx;k++){
				file<<std::fixed<<std::setprecision(3)<<U[i][j*nx+k]<<"\t"; 
			}
			file<<std::endl;
		}
		//file<<std::endl;
	}
	file.close();
	std::string ext2="_prop.txt";
	std::string fname2=filename+ext2;
	std::ofstream file2;
	file2.open(fname2);
	std::cout<<"Writing properties of wavefile file\n";
	file2<<std::fixed<<nx<<"\t"<<ny<<"\t"<<nt<<"\t"; 
	file2.close();
	
	std::cout<<"Writing is done\n";
	}

void write_rec_txt(const std::string& filename, double **U,int xsource,int ysource,double h,double dt , int nx,int ny, int nt){
	std::string ext=".txt";
	std::string fname=filename+ext;
	std::ofstream file;
	file.open(fname);
	
	for(int i=0;i<nt;i++){
		if(i%(nt/10)==0){
		std::cout<<"Writing Wavefield ... "<<float(i)/float(nt)*100<<"%\n";}
		for(int k=0;k<nx;k++){
				file<<std::fixed<<std::setprecision(3)<<U[i][k]<<"\t"; 
			}
		file<<std::endl;
	}
	file.close();
	std::string ext2="pr";
	std::ofstream file2;
	file2.open(filename+ext+ext2);
	file2<<nx<<"\t"<<nt<<"\t"<<h<<"\t"<<dt<<"\t"<<xsource<<"\t"<<ysource;
	file2<<std::endl;
	
	file2.close();
	std::cout<<"Writing is done\n";
	}

void write_img_txt(const std::string& filename, double **U,int nx,int ny){
	std::string ext=".txt";
	std::string fname=filename+ext;
	std::ofstream file;
	file.open(fname);
	
	std::cout<<"Writing Wavefield ... \n";	
	for(int i=0;i<ny;i++){
		for(int k=0;k<nx;k++){
				file<<std::fixed<<std::setprecision(3)<<U[i][k]<<"\t"; 
			}
		file<<std::endl;
	}
	file.close();
	std::cout<<"Writing is done\n";
	}

void read_rec(const std::string& filename, std::vector <double>& rec)
{
	std::ifstream f(filename);
	double content;
	while (f >> content)
	{
		rec.push_back(content);
	}
}

void read_rec_pr(const std::string& filename, int &nx,int &nt,int &xsource,int &ysource, double &h,double &dt)
{
	std::ifstream f(filename);
	int i=0;
	double content;
	while (f >> content)
	{
		i++;
		if(i==1){nx=content;}
		else if(i==2){nt=content;}
		else if(i==3){h=content;}
		else if(i==4){dt=content;}
		else if(i==5){ysource=content;}
		else if(i==6){xsource=content;}
	}
}

void rec_2d(double **rec_out,std::vector <double>& rec, int nx,int nt)
{
	for (int i=0;i<nt;i++)
	{
		for(int j=0;j<nx;j++)
		{
			rec_out[i][j]=rec[i*nx+j];
		}
	}
}

