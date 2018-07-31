import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.interpolate import interp2d 

vel=np.loadtxt("./VelocityModel/VelBP2004.50mx50m.txt")

nx=len(vel[0,:])
ny=len(vel[:,0])
L=nx*50.0
H=ny*50.0
x=np.linspace(0,L,nx)
y=np.linspace(0,H,ny)

f=interp2d(x,y,vel)
Lnew=L*1/3
xnew=np.arange(0,Lnew+20,20)
ynew=np.arange(0,H,20)
print len(xnew),len(ynew)
velnew=f(xnew,ynew)

traveltime=np.zeros([len(xnew)])

for i in range(len(xnew)):
	for j in range(len(ynew)):
		traveltime[i]+=20.0/velnew[j,i]

print "travel time max ",traveltime.max()

print "travel time min ",traveltime.min()
		
np.savetxt("./VelocityModel/VelBP2004.20mx20m.txt",velnew)
plt.figure(figsize=(8,3.5))
plt.title("Model BP2004",fontsize=12)
plt.imshow(velnew,extent=[0,Lnew/1000,H/1000,0])
plt.colorbar(orientation="horizontal").set_label("m/s")
plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
plt.savefig("Figures/vel_model.eps")
plt.show()
