
def fd2nd(u,h,nx,i):
    val=(-4*u[i]+u[i+1]+u[i-1]+u[i+nx]+u[i-nx])/h/h
    return val

def fd4th(u,h,nx,i):
    val=(-60*u[i]+16*u[i+1]+16*u[i-1]+16*u[i+nx]+16*u[i-nx]-u[i+2]-u[i-2]-u[i+2*nx]-u[i-2*nx])/h/h/12.0
    return val
 
def d2x(u,ux,h,nx,i):
	return (-2*u[i]+u[i+1]+u[i-1])*2.0/h/h-(ux[i+1]-ux[i-1])/2.0/h

def d2y(u,uy,h,nx,i):
	return (-2*u[i]+u[i+nx]+u[i-nx])*2.0/h/h-(uy[i+nx]-uy[i-nx])/2.0/h
	
def d3x(u,ux,h,nx,i):
	return (u[i+1]-u[i-1])*15.0/2.0/h/h/h-(ux[i+1]+8*ux[i]+ux[i-1])*3.0/2.0/h/h
	
def d3y(u,uy,h,nx,i):
	return (u[i+nx]-u[i-nx])*15.0/2.0/h/h/h-(uy[i+nx]+8*uy[i]+uy[i-nx])*3.0/2.0/h/h

def d2xy(u,ux,uy,h,nx,i):
	return (-2*uy[i]+uy[i-1]+uy[i+1])/h/h+(-ux[i-1-nx]-ux[i+1+nx]+ux[i-1]+ux[i+1]-2*(-2*ux[i]+ux[i-nx]+ux[i+nx]))/2.0/h/h+(5*u[i+1+nx]-5*u[i-1-nx]+u[i+1-nx]-u[i-1+nx]-4*u[i+nx]+4*u[i-nx]-6*u[i+1]+6*u[i-1])/4.0/h/h/h

def dx2y(u,ux,uy,h,nx,i):
	return (-2*ux[i]+ux[i-nx]+ux[i+nx])/h/h+(-uy[i-1-nx]-uy[i+1+nx]+uy[i-nx]+uy[i+nx]-2*(-2*uy[i]+uy[i-1]+uy[i+1]))/2.0/h/h+(5*u[i+1+nx]-5*u[i-1-nx]-u[i+1-nx]+u[i-1+nx]-4*u[i+1]+4*u[i-1]-6*u[i+nx]+6*u[i-nx])/4.0/h/h/h

def d4x(u,ux,h,nx,i):
    return -12.0/h/h/h/h*(u[i-1]-2*u[i]+u[i+1])+6.0/h/h/h*(ux[i+1]-ux[i-1])
    
def d4y(u,uy,h,nx,i):
    return -12.0/h/h/h/h*(u[i-nx]-2*u[i]+u[i+nx])+6.0/h/h/h*(uy[i+nx]-uy[i-nx])

def d2x2y(u,ux,uy,h,nx,i):
    return 1.0/h/h/h/h*(2.0*(u[i+nx]+u[i+1]-2.0*u[i]+u[i-1]+u[i-nx])- u[i-nx+1]-u[i-nx-1]-u[i+nx+1]-u[i+nx-1])+1.0/2.0/h/h/h*(ux[i+nx+1]+ux[i-nx+1]-ux[i-nx-1]-ux[i+nx-1]-2.0*ux[i+1]+2.0*ux[i-1])+1.0/2.0/h/h/h*(uy[i+nx+1]+uy[i+nx-1]-uy[i-nx-1]-uy[i-nx+1]-2.0*uy[i+nx]+2.0*uy[i-nx])
	
def d5x(u,ux,h,nx,i):
    return -90.0/h/h/h/h/h*(u[i+1]-u[i-1])+30.0/h/h/h/h*(ux[i+1]+4.0*ux[i]+ux[i-1])

def d5y(u,uy,h,nx,i):
    return -90.0/h/h/h/h/h*(u[i+nx]-u[i-nx])+30.0/h/h/h/h*(uy[i+nx]+4.0*uy[i]+uy[i-nx])
    
def d4xy(u,ux,h,nx,i):
    return -6.0/h/h/h/h/h*(u[i+nx+1]-u[i-1-nx]+u[i-1+nx]-u[i+1-nx]+2.0*u[i-nx]-2.0*u[i+nx])+3.0/h/h/h/h*(ux[i+nx+1]+ux[i-nx-1]-ux[i+nx-1]-ux[i-nx+1])

def dx4y(u,uy,h,nx,i):
    return -6.0/h/h/h/h/h*(u[i+nx+1]-u[i-1-nx]-u[i-1+nx]+u[i+1-nx]+2.0*u[i-1]-2.0*u[i+1])+3.0/h/h/h/h*(uy[i+nx+1]+uy[i-nx-1]-uy[i+nx-1]-uy[i-nx+1])
    
def d3x2y(u,ux,h,nx,i):
    return -3.0/2.0/h/h/h/h/h*(u[i+nx+1]-u[i-1-nx]+u[i+1-nx]-u[i-1+nx]+2.0*u[i-1]-2.0*u[i+1])+3.0/2.0/h/h/h/h*(ux[i+nx+1]+ux[i-nx-1]+ux[i+nx-1]+ux[i-nx+1]-2.0*ux[i+1]-2.0*ux[i-1])
	
def d2x3y(u,uy,h,nx,i):
    return -3.0/2.0/h/h/h/h/h*(u[i+nx+1]-u[i-1-nx]+u[i-1+nx]-u[i+1-nx]+2.0*u[i-nx]-2.0*u[i+nx])+3.0/2.0/h/h/h/h*(uy[i+nx+1]+uy[i-nx-1]+uy[i+nx-1]+uy[i-nx+1]-2.0*uy[i+nx]-2.0*uy[i-nx])
