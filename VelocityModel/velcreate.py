import numpy as np

vel1=np.ones([50,50])*4000
vel2=np.ones([100,100])*4000
vel3=np.ones([300,300])*4000
vel4=np.ones([500,500])*4000
vel5=np.ones([1000,1000])*4000
vel6=np.ones([3000,3000])*4000

np.savetxt("vel1.txt",vel1)
np.savetxt("vel2.txt",vel2)
np.savetxt("vel3.txt",vel3)
np.savetxt("vel4.txt",vel4)
np.savetxt("vel5.txt",vel5)
np.savetxt("vel6.txt",vel6)

