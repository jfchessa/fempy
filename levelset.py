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
    
    v=pt-n1
    
    A=np.array([a,b,n]) # A = [a;b;n]  transpose?
    xiv=np.linalg.solve(A.transpose(),v) # v/A;
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

    
def sdistq(pt,nodes):
    """
    Computes the signed distance to the midpoint of a triangle
    
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
    
    mp = np.average( nodes, axis=0 )
    
    v = pt-mp
    sdist = np.dot(v,n)
    
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
    for nid, pt in lspts.iteritems():
        dmin = 10.0e10
        
        for eid, elem in element.iteritems():
            tri = elem.Connectivity()
            de = sdistq(pt,node.CoordMat(tri)) #sdistpt2tri(pt,node.CoordMat(tri))
            
            if ( abs(de) < abs(dmin) ):
                dmin = de
                
        phi[n] = dmin
        n += 1
        
    return phi
        
#######################################################

# read in zero level set
zlsfile = io.GmshInput('fiber.msh')
zlsnodes = zlsfile.ReadNodes()
zlselem = zlsfile.ReadElements()


corners = zlsnodes.Corners()
corners[:,2] = 5*corners[:,2]

grid = msh.MeshHexa8(corners,11,11,7)
    
phi = mshgenls( zlsnodes, zlselem, grid.NodeArray() )

dataout = io.FeaData(grid.NodeArray(),grid.ElementArray())
dataout.AddNodeScalarField(phi,'Level Set')
dataout.WriteVtkFile('lstest')

dataout1 = io.FeaData(zlsnodes,zlselem)
dataout1.WriteVtkFile('zls')