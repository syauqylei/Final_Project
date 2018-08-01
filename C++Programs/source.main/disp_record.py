import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import sys

var=sys.argv[:]

rec1=np.loadtxt(var[1])
propadd=var[1]+"pr"
prop=np.loadtxt(propadd)
nx=int(prop[0])
nt=int(prop[1])
h=(prop[2])
dt=(prop[3])
x=np.linspace(0,nx*h,nx)
t=np.linspace(0,nt*dt,nt)
plt.pcolormesh(x,t,rec1[:][:],cmap="binary",vmax=100,vmin=-100)
plt.axis([0,x.max(),t.max(),0])
plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%d m'))
plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f s'))
plt.show()
