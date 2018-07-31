"""
This program show source function that is used in this research.
Based on Yang et al 2012.

"""

import matplotlib.pyplot as plt
import numpy as np
from  modules.source import *

t=np.arange(-0.25,0.25,0.001)
y=np.zeros(len(t))
f=30.0
for i in range(len(t)):
	y[i]=source_2(f,t[i])*0.0001

plt.plot(t,y)
plt.xlim(0,0.25)
plt.xlabel("Waktu (s)")
plt.ylabel("Normalisasi Amplitudo")
plt.savefig("Figures/sourcefunc.eps")
plt.show()
