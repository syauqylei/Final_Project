"""

this based on Yang et al 2012. "A central difference method with low dispersion effect for solving the scalar wave equation"

this program compare FD2nd , FD4th and NACD

"""

import matplotlib.pyplot as plt
import numpy as np
from modules.derv import *
from modules.wve import *

nx=301
ny=301
h=40
L=12000
H=12000
dt=0.005
fo=30
C=4000
nt=200

c=np.ones(nx*ny)*C

U2=wve_fd2nd([ny/2,nx/2],fo,c,nx,ny,nt,dt,h)
plt.imshow(U2[-1],cmap="binary",extent=[-6,6,-6,6])
plt.xticks(np.arange(-6,7,3))
plt.yticks(np.arange(-6,7,3))
plt.xlabel("Km")
plt.ylabel("Km")
plt.savefig("Figures/FD2nd.eps")
plt.show()
plt.clf()

U4=wve_fd4th([ny/2,nx/2],fo,c,nx,ny,nt,dt,h)
plt.imshow(U4[-1],cmap="binary",extent=[-6,6,-6,6])
plt.xlabel("Km")
plt.ylabel("Km")
plt.xticks(np.arange(-6,7,3))
plt.yticks(np.arange(-6,7,3))
plt.savefig("Figures/FD4th.eps")
plt.show()
plt.clf()

Unacd=wve_nacd([ny/2,nx/2],fo,c,nx,ny,nt,dt,h)
plt.imshow(Unacd[-1],cmap="binary",extent=[-6,6,-6,6])
plt.xlabel("Km")
plt.ylabel("Km")
plt.xticks(np.arange(-6,7,3))
plt.yticks(np.arange(-6,7,3))
plt.savefig("Figures/NACD.eps")
plt.show()
plt.clf()

plt.subplot(131)
plt.imshow(U2[-1],cmap="binary")
plt.axis('off')
plt.subplot(132)
plt.imshow(U4[-1],cmap="binary")
plt.axis('off')
plt.subplot(133)
plt.imshow(Unacd[-1],cmap="binary")
plt.axis('off')
plt.savefig("Figures/allmethods.eps")
plt.show()

