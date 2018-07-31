
def gen_in_sten_2nd(Nx,Ny):
    stencils=[]
    for i in range(Ny):
        for j in range(Nx):
            if i==0 or i==Ny-1 or j==0 or j==Nx-1:
                continue
            else:
                stencils.append(i*Nx+j)
    return stencils

def gen_in_sten_4th(Nx,Ny):
    stencils=[]
    for i in range(Ny):
        for j in range(Nx):
            if i==0 or i==1 or i==Ny-1 or i==Ny-2 or j==0 or j==1 or j==Nx-2 or j==Nx-1:
                continue
            else:
                stencils.append(i*Nx+j)
    return stencils

def gen_bdr_sten_2nd(Nx,Ny):
    stencils=[]
    for i in range(1,Ny-1):
        for j in range(1,Nx-1):
            if i==0 or i==Ny-1 or j==0 or j==Nx-1:
                stencils.append(i*Nx+j)
            else:
                continue
    return stencils
    
def gen_bdr_sten_4th(Nx,Ny):
    stencils=[]
    for i in range(Ny):
        for j in range(Nx):
            if i==0 or i==1 or i==Ny-1 or i==Ny-2 or j==0 or j==1 or j==Nx-2 or j==Nx-1:
                stencils.append(i*Nx+j)
            else:
                continue
    return stencils
