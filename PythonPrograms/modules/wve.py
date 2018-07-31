
import numpy as np
import matplotlib.pyplot as plt
from math import sin,exp, pi,sqrt
from derv import *
from stencils import *
from source import *

def wve_fd2nd(srcloc,f,c,nx,ny,nt,dt,h):
    print "Calculating Wavefield Using FD 2nd order ...."
    Nx=nx+2
    Ny=ny+2
    
    src_loc=srcloc[0]*Nx+1+srcloc[1]
    U=np.zeros([nt,ny,nx])
    u=np.zeros([nt+1,Nx*Ny])
    stens=gen_in_sten_2nd(Nx,Ny)
    for i in range(1,nt):
        t=dt*(i-1)
        for j in range(nx*ny):
                u[i+1,stens[j]]=2*u[i,stens[j]]-u[i-1,stens[j]]+(c[j]*dt)**2*(fd2nd(u[i],h,Nx,stens[j]))
        u[i+1,src_loc]+=source_2(f,t)
        
        for j in range(1,Ny-1):
            for k in range(1,Nx-1):
                U[i-1,j-1,k-1]=u[i,j*Nx+k]
        if i%(nt/10)==0:
            print i/float(nt)*100 ,"%"
	for i in range(1,Ny-1):
		for j in range(1,Nx-1):
			U[-1,i-1,j-1]=u[-1,i*Nx+j]
			
    return U

def wve_fd4th(srcloc,f,c,nx,ny,nt,dt,h):
    Nx=nx+4
    Ny=ny+4
    print "Calculating Wavefield Using FD 4th order ...."
    
    src_loc=srcloc[0]*Nx+2+srcloc[1]
    U=np.zeros([nt,ny,nx])
    u=np.zeros([nt+1,Nx*Ny])
    stens=gen_in_sten_4th(Nx,Ny)
    for i in range(1,nt):
        t=dt*(i-1)
        for j in range(nx*ny):
                u[i+1,stens[j]]=2*u[i,stens[j]]-u[i-1,stens[j]]+(c[j]*dt)**2*(fd4th(u[i],h,Nx,stens[j]))
        u[i+1,src_loc]+=source_2(f,t)
        for j in range(2,Ny-2):
            for k in range(2,Nx-2):
                U[i-1,j-2,k-2]=u[i,j*Nx+k]      
        if i%(nt/10)==0:
            print i/float(nt)*100 ,"%"
    for j in range(2,Ny-2):
		for k in range(2,Nx-2):
			U[-1,j-2,k-2]=u[-1,j*Nx+k]
    return U

def wve_nacd(srcloc,f,c,nx,ny,nt,dt,h):
    Nx=nx+2
    Ny=ny+2
    print "Calculating Wavefield Using Nearly Analytic Discrete Method ...."
    
    src_loc=srcloc[0]*Nx+1+srcloc[1]
    U=np.zeros([nt,ny,nx])
    u=np.zeros([nt+1,Nx*Ny])
    ux=np.zeros([nt+1,Nx*Ny])
    uy=np.zeros([nt+1,Nx*Ny])
    stens=gen_in_sten_2nd(Nx,Ny)
    
    for i in range(1,nt):
        t=dt*i
        for j in range(nx*ny):
                u[i+1,stens[j]]=2*u[i,stens[j]]-u[i-1,stens[j]]+(c[j]*dt)**2*(fd2nd(u[i],h,Nx,stens[j]))-((c[j]*h*dt)**2-((c[j]*dt)**4))/12.0*(d4x(u[i],ux[i],h,Nx,stens[j])+d4y(u[i],uy[i],h,Nx,stens[j]))+((c[j]*dt)**4)/6.0*(d2x2y(u[i],ux[i],uy[i],h,Nx,stens[j]))
                ux[i+1,stens[j]]=2*ux[i,stens[j]]-ux[i-1,stens[j]]+(c[j]*dt)**2*(fd2nd(ux[i],h,Nx,stens[j]))-((c[j]*h*dt)**2-((c[j]*dt)**4))/12.0*(d5x(u[i],ux[i],h,Nx,stens[j])+dx4y(u[i],uy[i],h,Nx,stens[j]))+((c[j]*dt)**4)/6.0*(d3x2y(u[i],ux[i],h,Nx,stens[j]))
                uy[i+1,stens[j]]=2*uy[i,stens[j]]-uy[i-1,stens[j]]+(c[j]*dt)**2*(fd2nd(uy[i],h,Nx,stens[j]))-((c[j]*h*dt)**2-((c[j]*dt)**4))/12.0*(d4xy(u[i],ux[i],h,Nx,stens[j])+d5y(u[i],uy[i],h,Nx,stens[j]))+((c[j]*dt)**4)/6.0*(d2x3y(u[i],uy[i],h,Nx,stens[j]))
        u[i+1,src_loc]+=source_2(f,t)
        for j in range(ny):
            for k in range(nx):
                U[i-1,j,k]=u[i,stens[j*nx+k]]
        if i%(nt/10)==0:
            print i/float(nt)*100 ,"%"

    
	for j in range(ny):
		for k in range(nx):
			U[-1,j,k]=u[-1,stens[j*nx+k]]
    return U
