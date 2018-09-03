import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import sys

var=sys.argv[:]

rec1=np.loadtxt(var[1])
nx=450
ny=175
h=20
dt=20
x=np.linspace(0,nx*h,nx)
t=np.linspace(0,ny*dt,ny)
plt.pcolormesh(x,t,rec1[:][:],cmap="binary")
plt.axis([0,x.max(),t.max(),0])
plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%d m'))
plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%d m'))
plt.show()
