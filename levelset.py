import numpy as np
import fempy as fp
import fempy.meshing as msh
import fempy.io as io

def distpt2line(pt,n1,n2):
    """
     [dist,cp]=distpt2line(pt,n1,n2)

    Computes the distance, d, from a point, pt to a finite line segment 
    defined by its endpoints given by n1 and n2, as well as the point on the 
    line segment that is the closest to the point, cp.  These points must be 
    given by a row vector of length 3

    written by: Jack Chessa, jfchessa@utep.edu
    """
    
    u=n2-n1
    v=pt-n1
    
    udv = np.dot(u,v)
    lv = np.linalg.norm(v)
    
    if ( udv<=0 ):  # cp=n1
        d=lv
        cp=n1
        
    else:
    
        lu = np.linalg.norm(u)
        lu2 = lu**2
    
        if ( udv>=lu2 ):   # cp = n2 
            cp=n2
    
        else:  # cp is on the ine segment
            cp=n1+(udv/lu2)*u
        
        d=np.linalg.norm(pt-cp)
        
    return [d,cp]
    
def sdistpt2tri(pt,nodes):
    """
    Computes the signed distance and closest point from a point to a finite
    triangle.
    
        SDIST -  the signed distance
        CP -  the coordinates of the closest point on the finite triangle to
              the point
        PT - the row vector containing the coordiants of the point
        N1, N2, N3 - the row vectors containing the coordiants of the verticies 
                     of the triangle.  
  
    The sign convention is such that the sign is positive if the point is
    in the positive direction w.r.t to the outward normal of the triangle.
    The outward normal is defined in the right handed sense of the node
    numbering of N1, N2, N3.

    written by: Jack Chessa, jfchessa@utep.edu
    """
    n1 = nodes[0]
    n2 = nodes[1]
    n3 = nodes[2]
    
    a=n2-n1
    b=n3-n1
    n=np.cross(a,b)
    n=n/np.linalg.norm(n)
    #n4=n1+n
    
    v=pt-n1
    
    A=np.array([a,b,n]) # A = [a;b;n]
    xiv=np.linalg.solve(A,v) # v/A;
    A1=xiv[0]
    A2=xiv[1]
    A3=1-A1-A2
    zeta=xiv[2]
    
    if ( zeta<0 ):
        sgn=-1
    else:
        sgn=1
    
    if ( A1<=0 ):
        if ( A2<=0 ):  # cp is n1
            cp=n1
            
        elif ( A3<=0 ): # n3 is cp
            cp=n3
            
        else:  # cp is on b
            [d,cp]=distpt2line(pt,n1,n3)
            sdist=d*sgn
            return sdist
            
        sdist=np.linalg.norm(pt-cp)*sgn
        
    elif ( A2<=0 ):
        if ( A3<=0 ): # cp is n2
            cp=n2
            
        else: # cp is on a
            [d,cp]=distpt2line(pt,n1,n2)
            sdist=d*sgn
            return sdist
            
        sdist=np.linalg.norm(pt-cp)*sgn
        
    elif (A3<0 ): # cp is on c
        [d,cp]=distpt2line(pt,n3,n2)
        sdist=d*sgn
        return sdist
        
    else: # point over triangle
        sdist=zeta
        dvect=sdist*n
        cp=pt-dvect
        
    return sdist

    
def mshgenls(node,element,lspts):
    """
    Generates a level set function from a triagulated mesh of a zero
    levelset.

      NODE - a node coordinate matrix for the zero level set mesh
      ELEMENT - an element connectivity matrix for the zero level set mesh
                (must be tri3)
      LSPTS - a node coordinate matrix of the points where the level set is
              to be computed
    """
    phi = np.zeros( len(lspts), dtype=float )
    n=0
    for pt in lspts:
        dmin = 10.0e10
        
        for tri in element:
            de = sdistpt2tri(pt,node[tri])
            
            if ( abs(de) < abs(dmin) ):
                dmin = de
                
        phi[n] = dmin
        n += 1
        
    return phi
        
#######################################################
L=10
W=10
H=10

corners = np.array([[0,0,0],[L,0,0],[L,W,0],[0,W,0],
                    [0,0,H],[L,0,H],[L,W,H],[0,W,H]])
nodes = msh.node_array3d(corners,2,2,2)

zlscorners = np.array([[0,0,H/2],[L,0,H/2],[L,W,H/2],[0,W,H/2]])
nnx=2
nny=2
zlsnodes = msh.node_array2d(zlscorners,nnx,nny)
zlsconn = np.concatenate ( ( msh.gen_conn2d(np.array([0,1,nnx],dtype=int),nnx-1,nny-1),
      msh.gen_conn2d(np.array([1,nnx+1,nnx],dtype=int),nnx-1,nny-1) ), axis=0 )
      
phi = mshgenls( zlsnodes, zlsconn, nodes )

dataout = io.FeaData(nodes)
dataout.SetDisplacement(dd,dofmap)
dataout.SetStress(stress)
dataout.AddCellScalarField(svm,"mises")
dataout.WriteVtkFile(filename)