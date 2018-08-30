import numpy as np

a=np.zeros([175,450])
txt=".txt"
for i in range(20,440,20):
	print i
	b=np.genfromtxt("image"+str(i)+txt)
	a+=b

np.savetxt("image.txt",a)
