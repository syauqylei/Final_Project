#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

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

void write_rec_txt(const std::string& filename, double **U,int xsource,int ysource,double h,double dt , int nx, int nt){
	std::string fname=filename;
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
	file2.open(filename+ext2);
	file2<<nx<<"\t"<<nt<<"\t"<<h<<"\t"<<dt<<"\t"<<ysource<<"\t"<<xsource;
	file2<<std::endl;
	
	file2.close();
	std::cout<<"Writing is done\n";
	}
