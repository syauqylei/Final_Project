This repository is used for my final project research. The research focus on simulation of Wavefield used on seismic proccessing.
Basic Idea
With NACD method proposed by Yang et al 2012 created a way to simulate wavefield modelling with less dispersion error that could make
seismic proccessing method like RTM become more accurate to implement

PythonPrograms/
this Directory consist of programs that compare NACD wavefield simulation among of Finite difference simulation with different order of accurations.
it uses python 2.7 version.

Error calculation this research use same way that is used by Yang et al, 2006.

to see Errors Comparison

$python PythonPrograms/error_cal.py

to see Source Function

$python PythonPrograms/Source_function.py 

to see Visual comparison # Take long time to run about 30 mnt in my i5 laptop with 4gb ram

$python PythonPrograms/2d_Vis_Compar.py

to see Velocity model that is used (Velocity Model B2004)

$python PythonPrograms/Vel_read.py  

all figures will be saved directly to Figures/ Directory
