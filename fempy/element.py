import numpy as np
import math

import basic
import quadrature as quad

import pdb

def sctr_array(econn,ndofn):
    nn = len(econn)
    ss = nn*ndofn
    sctr = np.ndarray( ss, basic.INDX_TYPE )
    s0 = ndofn*np.array(econn,basic.INDX_TYPE)
    for s in xrange(ndofn):
        sctr[s:ss:ndofn] = s0 + s
    return sctr
    
def elem_jacobian(dNdxi,coord,sdim=None,edim=None):
    """computes the element jacobian matrix for element"""
    
    if ( coord.ndim == 1 ): 
        coord = np.reshape( coord, (coord.size,1) )
    
    if ( sdim == None ):
        sdim = coord.shape[1]
        
    if ( edim == None ):
        edim = dNdxi.shape[1]
    
    jmat = np.zeros( (sdim,sdim), dtype=basic.FLOAT_TYPE )
    
    if ( edim == 3 ):
        jmat[:sdim,:edim] = np.dot( coord[:,:sdim].T, dNdxi[:,:edim] )
        return jmat
        
    if ( edim == 2 ):
        jmat[:sdim,:edim] = np.dot( coord[:,:sdim].T, dNdxi[:,:edim] )
        if ( sdim == 3 ):  
            jmat[:,2] = np.cross( jmat[:,0], jmat[:,1] )
            jmat[:,2] = jmat[:,2]/np.linalg.norm(jmat[:,2])
        
        return jmat
            
    else: #  edim = 1 
        if ( sdim == 3 ):
            jmat[:sdim,:edim] = np.dot( coord[:,:sdim].T, dNdxi[:,:edim] )
            jmat[:,1] = np.cross( np.array([0,0,1],dtype=basic.FLOAT_TYPE), jmat[:,0] )
            if ( np.linalg.norm(jmat[:,1])==0 ):
                jmat[:,1] = np.cross( np.array([0,1,0],dtype=basic.FLOAT_TYPE), jmat[:,0] )
            jmat[:,1] = jmat[:,1]/np.linalg.norm(jmat[:,1])
            jmat[:,2] = np.cross( jmat[:,0], jmat[:,1] )
            jmat[:,2] = jmat[:,2]/np.linalg.norm(jmat[:,2])
            
        elif ( sdim == 2 ):
            jmat[:sdim,:edim] = np.dot( coord[:,:sdim].T, dNdxi[:,:edim] )
            jmat[0,1] = -jmat[1,0]
            jmat[1,1] = jmat[0,0]
            jmat[:,1] = jmat[:,1]/np.linalg.norm(jmat[:,1])
            
        else: # sdim = 1
            jmat = np.dot( coord[:,:sdim].T, dNdxi )
        
        return jmat
    
def grad_basis(dNdxi,coord,sdim=None,edim=None):
    
    if ( coord.ndim == 1 ): 
        coord = np.reshape( coord, (coord.size,1) )
        
    if ( sdim == None ):
        sdim = coord.shape[1]
        
    if ( edim == None ):
        edim = dNdxi.shape[1]
        
    jac = elem_jacobian(dNdxi,coord,sdim,edim)
    djac = np.linalg.det(jac)  # there is a bit of a redundancy here
    invj = np.linalg.inv(jac)
    
    if ( not(edim==1) ):
        dNdx = np.dot(dNdxi,invj[:edim,:])
        return [dNdx,djac]
    dNdx = np.dot(dNdxi,invj[:edim,:])
    
    return [dNdx,djac]
        
        

def form_bmat_2d(dndx):
    nn=dndx.shape[0]
    bmat=np.zeros([3,2*nn], dtype=basic.FLOAT_TYPE)
    bmat[0,0:2*nn:2] = dndx.transpose()[0,:]
    bmat[1,1:2*nn:2] = dndx.transpose()[1,:]
    bmat[2,0:2*nn:2] = dndx.transpose()[1,:]
    bmat[2,1:2*nn:2] = dndx.transpose()[0,:]
    return bmat
    
    
def form_bmat_3d(dndx):
    nn=dndx.shape[0]
    bmat=np.zeros([6,3*nn], dtype=basic.FLOAT_TYPE)
    bmat[0,0:3*nn:3] = dndx.transpose()[0,:]
    bmat[1,1:3*nn:3] = dndx.transpose()[1,:]
    bmat[2,2:3*nn:3] = dndx.transpose()[2,:]
    bmat[5,0:3*nn:3] = dndx.transpose()[1,:]
    bmat[5,1:3*nn:3] = dndx.transpose()[0,:]
    bmat[3,1:3*nn:3] = dndx.transpose()[2,:]
    bmat[3,2:3*nn:3] = dndx.transpose()[1,:]
    bmat[4,2:3*nn:3] = dndx.transpose()[0,:]
    bmat[4,0:3*nn:3] = dndx.transpose()[2,:]
    return bmat

#######################################################################
#
#  1 D   E L E M E N T S
#

class ElemLine2(object):
    
    def __init__(self,carray=[],prop=[]):
        self.prop = prop
        assert ( len(carray) == self.NumNodes() ), 'Incorrect number of nodes in connectivity'
        self.conn = np.array(carray,dtype= INDX_TYPE)
        
    def __repr__(self):
        return 'Line2 element '+str(self.conn)
        
    def Type(self):
        return 'LINE2'
        
    def NumNodes(self):
        return 2

    def Connectivity(self):
        return self.conn
        
    def NodeCoord(self,nodes):
        return nodes[self.conn]
        
    def Property(self):
        return self.prop
        
    def Material(self):
        return self.prop['Matl']
        
    def ElemDim(self):
        return 1
        
    def Order(self):
        return 1
        
    def Topology(self):
        return 'Line'
        
    def QuadratureRule(self,p=1):
        n = int( math.ceil(0.5*(p+1)))
        qpt, qwt =  quad.quadrature_gauss1d(n)
        return [ np.reshape(qpt,(len(qpt),1)), qwt ]
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype=basic.FLOAT_TYPE )
        basis = (0.5)*np.array( [(1-pt[0]), (1+pt[0]) ], dtype=basic.FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype=basic.FLOAT_TYPE )
        dbasis = np.array([[ -0.5 ],[ 0.5 ]], dtype=basic.FLOAT_TYPE )
        return dbasis

    def dNdx(self,coords,pt=None,sdim=None):
        
        if ( coords.ndim == 1 ): 
            coords = np.reshape( coords, (coords.size,1) )
        
        if ( sdim == None ):
            sdim = coords.shape[1]
        
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype=basic.FLOAT_TYPE )
              
        return grad_basis(self.dNdxi(pt),coords,sdim=sdim,edim=self.ElemDim())
      
    def Jacobian(self,coords,pt=None,sdim=None):
        jmat = elem_jacobian(self.dNdxi(pt),coords,sdim=sdim,edim=self.ElemDim())
        return np.linalg.det(jmat)
              
    def BMat(self,coords,pt=None,sdim=None):
        
        if ( coords.ndim == 1 ): 
            coords = np.reshape( coords, (coords.size,1) )
        
        if ( sdim == None ):
            sdim = coords.shape[1]
        
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype=basic.FLOAT_TYPE )
            
        dNdx = self.dNdx(coords,pt,sdim)
        
        if ( self.ElemDim() == sdim ):
            return dNdx
            
        bmat = np.zeros( (1,sdim*self.NumNodes()), dtype=basic.FLOAT_TYPE )
        for s in xrange( sdim ):
            bmat[0,s:bmat.size:sdim] = dNdx[0][:,s]
        
        return [ bmat, dNdx[1] ]   

class ElemTruss3D(ElemLine2):
     
    #def __init__(self,carray=[],prop={}):
    #    self.prop = prop
    #    assert ( len(carray) == self.NumNodes() ), 'Incorrect number of nodes in connectivity'
    #    self.conn = np.array(carray,dtype= INDX_TYPE)
        
    def __repr__(self):
        return '3D truss element '+str(self.conn)
        
    def Type(self):
        return 'TRUSS3D'


#######################################################################
#
#  2 D   E L E M E N T S
#
    
class ElemTria3(ElemLine2):
    
    #def __init__(self,carray,prop=[]):
    #    self.prop = prop
    #    assert ( len(carray) == self.NumNodes() ), 'Incorrect number of nodes in connectivity'
    #    self.conn = np.array(carray,dtype= INDX_TYPE)
        
    def __repr__(self):
        return 'Tria3 element '+str(self.conn)
        
    def Type(self):
        return 'TRIA3'
        
    def NumNodes(self):
        return 3
        
    def ElemDim(self):
        return 2
        
    def Order(self):
        return 1
        
    def Topology(self):
        return 'Tria'
        
    def QuadratureRule(self,p=1):
        return  quad.quadrature_simplex(self.ElemDim(),p)
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        basis = np.array( [1-pt[0]-pt[1], pt[0], pt[1] ], dtype= basic.FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        dbasis = np.array([ [-1,-1], [1,0], [0,1] ], dtype= basic.FLOAT_TYPE )
        return dbasis
        
    def BMat(self,coords,pt=None,sdim=None):
        
        if ( coords.ndim == 1 ): 
            coords = np.reshape( coords, (coords.size,1) )
        
        if ( sdim == None ):
            sdim = coords.shape[1]
        
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
            
        ipm = self.dNdx(coords,pt,sdim)
        
        if ( self.ElemDim() == sdim ):
            return [ form_bmat_2d(ipm[0]), ipm[1] ]     
 
        bmat = np.zeros( (6,sdim*self.NumNodes()), dtype= basic.FLOAT_TYPE )
        bmat[0,0:bmat.size:sdim] = ipm[0][:,0]
        bmat[1,1:bmat.size:sdim] = ipm[0][:,1]
        bmat[2,2:bmat.size:sdim] = ipm[0][:,2]
        
        bmat[5,1:bmat.size:sdim] = ipm[0][:,2]
        bmat[5,2:bmat.size:sdim] = ipm[0][:,1]
        
        bmat[3,0:bmat.size:sdim] = ipm[0][:,2]
        bmat[3,2:bmat.size:sdim] = ipm[0][:,0]
        
        bmat[4,0:bmat.size:sdim] = ipm[0][:,1]
        bmat[4,1:bmat.size:sdim] = ipm[0][:,0]
        
        return [ bmat, ipm[1] ]   
    
class ElemQuad4(ElemTria3):
    
    #def __init__(self,carray,prop=[]):
    #    self.prop = prop
    #    self.conn = np.array(carray,dtype= INDX_TYPE)
        
    def __repr__(self):
        return 'Quad4 element '+str(self.conn)
        
    def Type(self):
        return 'QUAD4'
        
    def NumNodes(self):
        return 4
        
    def Order(self):
        return 1
        
    def Topology(self):
        return 'Quad'
        
    def QuadratureRule(self,p=1):
        n = int( math.ceil(0.5*(p+1)))
        qr =  quad.quadrature_gauss1d(n)
        return  quad.compound_quadrature(qr[0],qr[1],qr[0],qr[1])
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        basis = (0.25)*np.array( [(1-pt[0])*(1-pt[1]), \
                     (1+pt[0])*(1-pt[1]), \
                     (1+pt[0])*(1+pt[1]), \
                     (1-pt[0])*(1+pt[1]) ], dtype= basic.FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        dbasis = (0.25)*np.array([ [ -(1-pt[1]), -(1-pt[0]) ],\
                                   [  (1-pt[1]), -(1+pt[0]) ],\
                                   [  (1+pt[1]),  (1+pt[0]) ],\
                                   [ -(1+pt[1]),  (1-pt[0]) ]], dtype= basic.FLOAT_TYPE )
        return dbasis
   
 

#######################################################################
#
#  3 D   E L E M E N T S
#

class ElemTetra4(ElemTria3):
    
    #def __init__(self,carray,prop=[]):
    #    self.prop = prop
    #    self.conn = np.array(carray,dtype= INDX_TYPE)
        
    def __repr__(self):
        return 'Tetra4 element '+str(self.conn)
        
    def Type(self):
        return 'TETRA4'
        
    def NumNodes(self):
        return 4
        
    def ElemDim(self):
        return 3
        
    def Order(self):
        return 1
        
    def Topology(self):
        return 'Tetra'
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        basis = np.array( [1-pt[0]-pt[1]-pt[2], pt[0], pt[1], pt[2] ], dtype= basic.FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        dbasis = np.array([ [-1,-1,-1], [1,0,0], [0,1,0], [0,0,1] ], dtype= basic.FLOAT_TYPE )
        return dbasis
        
    def BMat(self,coords,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        ipm = self.dNdx(coords,pt)
        return [ form_bmat_3d(ipm[0]), ipm[1] ]     
 
#----------- HEXA8 Element ------------
class ElemHexa8(ElemTetra4):
    
    #def __init__(self,carray,prop=[]):
    #    self.prop = prop
    #    self.conn = np.array(carray,dtype= INDX_TYPE)
        
    def __repr__(self):
        return 'Hexa8 element '+str(self.conn)
        
    def Type(self):
        return 'HEXA8'
        
    def NumNodes(self):
        return 8
        
    def Order(self):
        return 1
        
    def Topology(self):
        return 'Hexa'
        
    def QuadratureRule(self,p=1):
        n = int( math.ceil(0.5*(p+1)))
        qr =  quad.quadrature_gauss1d(n)
        return  quad.compound_quadrature(qr[0],qr[1],qr[0],qr[1],qr[0],qr[1])
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        basis = .125*np.array( [(1-pt[0])*(1-pt[1])*(1-pt[2]), \
                     (1+pt[0])*(1-pt[1])*(1-pt[2]), \
                     (1+pt[0])*(1+pt[1])*(1-pt[2]), \
                     (1-pt[0])*(1+pt[1])*(1-pt[2]), \
                     (1-pt[0])*(1-pt[1])*(1+pt[2]), \
                     (1+pt[0])*(1-pt[1])*(1+pt[2]), \
                     (1+pt[0])*(1+pt[1])*(1+pt[2]), \
                     (1-pt[0])*(1+pt[1])*(1+pt[2]) ], dtype= basic.FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= basic.FLOAT_TYPE )
        dbasis = 0.125*np.array([ [-(1-pt[1])*(1-pt[2]), -(1-pt[0])*(1-pt[2]), -(1-pt[0])*(1-pt[1]) ],\
                                   [  (1-pt[1])*(1-pt[2]), -(1+pt[0])*(1-pt[2]), -(1+pt[0])*(1-pt[1]) ],\
                                   [  (1+pt[1])*(1-pt[2]),  (1+pt[0])*(1-pt[2]), -(1+pt[0])*(1+pt[1]) ],\
                                   [ -(1+pt[1])*(1-pt[2]),  (1-pt[0])*(1-pt[2]), -(1-pt[0])*(1+pt[1]) ],\
                                   [ -(1-pt[1])*(1+pt[2]), -(1-pt[0])*(1+pt[2]),  (1-pt[0])*(1-pt[1]) ],\
                                   [  (1-pt[1])*(1+pt[2]), -(1+pt[0])*(1+pt[2]),  (1+pt[0])*(1-pt[1]) ],\
                                   [  (1+pt[1])*(1+pt[2]),  (1+pt[0])*(1+pt[2]),  (1+pt[0])*(1+pt[1]) ],\
                                   [ -(1+pt[1])*(1+pt[2]),  (1-pt[0])*(1+pt[2]),  (1-pt[0])*(1+pt[1]) ] ], dtype= basic.FLOAT_TYPE )
        return dbasis
 
    	
class ElementArray(dict):
    
    def __init__(self,*args):
        self.elems = {}

 
    def __repr__(self):
        s='Element Array\n'
        for eid, e in self.iteritems():
            s=s+str(eid)+': ' +str(e)+'\n'
        return s
        
###############################################################################        
        
class ElementArray(dict):
    
    def __init__(self,*args):
        self.elems = {}

 
    def __repr__(self):
        s='Element Array\n'
        for eid, e in self.iteritems():
            s=s+str(eid)+': ' +str(e)+'\n'
        return s
        
               	    	    	              