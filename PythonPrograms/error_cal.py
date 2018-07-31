#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 12:14:06 2018

@author: syauqy
this program based on Yang et al 2006. "Optmal Nearly Analytic Discrete Approximation to the scalar wave equation"
"""
#Calculation of errors
import numpy as np
import matplotlib.pyplot as plt
from modules.derv import *
from math import cos,sin,sqrt,pi

h=25
nt=2000
nx=205
fo=20
dt=0.001
x=np.linspace(-50,5050,nx)
U=np.zeros([nt,205])
c=4000
t=0
U=np.zeros([nt,nx])
Ux=np.zeros([nt,nx])
U2=np.zeros([nt+1,nx])
U4=np.zeros([nt+1,nx])
Unacd=np.zeros([nt+1,nx])
Unacdx=np.zeros([nt+1,nx])

U0=np.zeros([nx])
U0x=np.zeros([nx])
Ut=np.zeros([nx])
Utx=np.zeros([nx])
Um1=np.zeros([nx])
Um1x=np.zeros([nx])

U=np.zeros([nt,nx])

#Define Initial condition
for i in range(nx):
    U0[i]=cos(-2*pi*fo*x[i]/c)
    U0x[i]=2*pi*fo/c*sin(-2*pi*fo*x[i]/c)
for i in range(nx):
    Ut[i]=-2*pi*fo*sin(-2*pi*fo*x[i]/c)
    Utx[i]=(2*pi*fo)**2/c*cos(-2*pi*fo*x[i]/c)


Um1=U0-Ut*dt
Um1x=U0x-Utx*dt
for i in range(nt):
    t=dt*i
    for j in range(nx):
        U[i,j]=cos(2*pi*fo*(t-x[j]/c)) 
        Ux[i,j]=2*pi*fo/c*sin(2*pi*fo*(t-x[j]/c))
t=0
U2[0]=Um1
U2[1]=U0
U4[0]=Um1
U4[1]=U0
Unacd[0]=Um1
Unacd[1]=U0
Unacdx[0]=Um1x
Unacdx[1]=U0x

#calculting Wavefield according to dinghui yang 2012 NACD
for i in range(1,nt-1):
    t=dt*(i-1)
    U2[i,0]=U[i-1,0]
    U2[i,-1]=U[i-1,-1]
    Unacd[i,0]=U[i-1,0]
    Unacd[i,-1]=U[i-1,-1]
    Unacdx[i,0]=Ux[i-1,0]
    Unacdx[i,-1]=Ux[i-1,-1]
    for j in range(1,nx-1):
        U2[i+1,j]=2*U2[i,j]-U2[i-1,j]+(c*dt/h)**2*(U2[i,j+1]+U2[i,j-1]-2*U2[i,j])
        Unacd[i+1,j]=2*Unacd[i,j]-Unacd[i-1,j]+(c*dt/h)**2*(Unacd[i,j+1]+Unacd[i,j-1]-2*Unacd[i,j])-((c*h*dt)**2-(c*dt)**4)/12.0*(d4x(Unacd[i],Unacdx[i],h,nx,j))
        Unacdx[i+1,j]=2*Unacdx[i,j]-Unacdx[i-1,j]+(c*dt/h)**2*(Unacdx[i,j+1]+Unacdx[i,j-1]-2*Unacdx[i,j])-((c*h*dt)**2-(c*dt)**4)/12.0*(d5x(Unacd[i],Unacdx[i],h,nx,j))
    U4[i,:2]=U[i-1,:2]
    U4[i,-2:]=U[i-1,-2:]
    for j in range(2,nx-2):
        U4[i+1,j]=2*U4[i,j]-U4[i-1,j]+(c*dt/h)**2*(16*U4[i,j+1]+16*U4[i,j-1]-U4[i,j-2]-U4[i,j+2]-30*U4[i,j])/12
U2[-2,0]=U[-1,0]
U2[-2,-1]=U[-1,-1]
Unacd[-2,0]=U[-1,0]
Unacd[-2,-1]=U[-1,-1]
U4[-2,:2]=U[-1,:2]
U4[-2,-2:]=U[-1,-2:]

err1=np.zeros(nt-1)
err2=np.zeros(nt-1)
err3=np.zeros(nt-1)
for i in range(nt-1):
    err1[i]=sqrt(sum((U[i]-U2[i+1])**2)/sum(U[i]**2))*100
    err2[i]=sqrt(sum((U[i]-Unacd[i+1])**2)/sum(U[i]**2))*100
    err3[i]=sqrt(sum((U[i]-U4[i+1])**2)/sum(U[i]**2))*100

xx=np.linspace(0,2,nt-1)
plt.semilogy(xx,err1,"r-",label="FD2",linewidth=0.4)
plt.semilogy(xx,err3,"g-",label="FD4",linewidth=0.4)
plt.semilogy(xx,err2,"b-",label="NACD",linewidth=0.4)
plt.legend()
plt.ylim(1e-1,1e+2)
plt.xlim(0,2)
plt.xlabel("Waktu (s)")
plt.ylabel("Error Relatif (%)")
plt.savefig("Figures/Graphic_of_errors.eps")
plt.show()
