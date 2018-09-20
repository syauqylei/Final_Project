import numpy as np
import matplotlib.pyplot as plt

cpu=np.array([0.178,0.872,5.438,14.133,57.304])
gpu=np.array([0.132,0.151,0.541,1.152,4.324])
gpum=np.array([0.993,1.766,0.641,1.246,4.498])

elms=np.array([50*50,100*100,300*300,500*500,1000*1000])

plt.semilogx(elms,cpu,"r-",label="CPU")
plt.semilogx(elms,gpu,"g-",label="GPU")
plt.semilogx(elms,gpum,"b-",label="GPU-managed")
plt.xlim(2500,1e+6)
plt.legend()
plt.xlabel("Jumlah elemen")
plt.ylabel("Waktu (s)")
plt.show()
