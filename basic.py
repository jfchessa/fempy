import numpy as np

#import scipy as sp
from scipy import sparse
from scipy.sparse.linalg import dsolve

import pdb

# numerical type definitions
FLOAT_TYPE  =  np.float64   # floats
#INDX_TYPE   =  np.uint32    # unsigned int
INDX_TYPE   =  np.int32     # ints for index things
        
        
def scatter_matrix(ke,Kmat,sctr1,sctr2=None):
	if ( sctr2==None ):
		sctr2=sctr1
	for i in xrange(len(sctr1)):
		ii=sctr1[i]
		for j in xrange(len(sctr2)):
			jj=sctr2[j]
			Kmat[ii,jj] += ke[i,j]

class DelayedAssm(object):
    
    def __init__(self,p1,p2=None): 
        
        if ( p2==None ):
            expdim = p1
        else:
            expdim = self.CalcSpace(p1,p2)
              
        self.Kg = np.zeros( expdim, FLOAT_TYPE )
        self.Ig = np.zeros( expdim, INDX_TYPE )
        self.Jg = np.zeros( expdim, INDX_TYPE )
        
        ndfeMax = 60
        self.ii = np.outer(range(ndfeMax),np.ones(ndfeMax,INDX_TYPE))
        self.jj = self.ii.T
        self.kk = np.array( range(ndfeMax**2), INDX_TYPE )
      
    def CalcSpace(self,element,dpn):
        expdim = 0
        for e in element:
            expdim = expdim + (e.NumNodes()*dpn)**2
        return expdim 
         
    def AddLocalMatrix(self,ke,sctr,sctr2=None):
        
        if sctr2==None:  # square matrix
            ndfe = len(sctr)
            ndfe2 = ndfe**2
            self.Ig[self.kk[:ndfe2]] = sctr[self.ii[:ndfe,:ndfe]]
            self.Jg[self.kk[:ndfe2]] = sctr[self.jj[:ndfe,:ndfe]]
            self.Kg[self.kk[:ndfe2]] = ke[:ndfe,:ndfe]
            self.kk = self.kk + ndfe2
            
        else:  # none square matrix
            ndfei = len(sctr)
            ndfej = len(sctr2)
            ndfe2 = ndfei*ndfej
            self.Ig[self.kk[:ndfe2]] = sctr[self.ii[:ndfei,:ndfej]]
            self.Jg[self.kk[:ndfe2]] = sctr2[self.jj[:ndfei,:ndfej]]
            self.Kg[self.kk[:ndfe2]] = ke[:ndfei,:ndfej]
            self.kk = self.kk + ndfe2
  
    def GetCsrMatrix(self):
    
        return sparse.csr_matrix( (self.Kg, (self.Ig,self.Jg)) )
        self.Kg = 0
        self.Ig = 0
        self.Jg = 0
       
    def NumDof(self):
        return max(self.Ig)
         
        
        
        
# this seems to be quite a bit slower than fesolve.  It zeros out the
# rows and columns in the matrix
def fesolve2(Kmat,fext,ispc,vspc=None):
	n=np.shape(Kmat)[0]
	tr=0.0
	for i in xrange(n):
		tr += Kmat[i,i]
	beta=tr/(1.0*n)
	nspc=len(ispc)
	
	Kreac=Kmat[ispc,:].tocsr()
	if vspc==None:
		vspc=np.zeros((nspc,1))
	else: # inhomogenious bc
		for i in xrange(nspc):
			ii = ispc[i]
			iv = vspc[i]
			for jj in xrange(Kmat.shape[0]):
			     fext[jj] = fext[jj] - Kmat[jj,ii]*iv
			
	for i in xrange(nspc):
  	       ii = ispc[i]
  	       for jj in xrange(n):
  	           if ( not(Kmat[ii,jj]==0) ):
	               Kmat[ii,jj] = 0.0
	           if ( not(Kmat[jj,ii]==0) ):
	               Kmat[jj,ii] = 0.0
	       Kmat[ii,ii] = beta
	       fext[ii] = beta*vspc[i]

	Kmat=Kmat.tocsr()
	d = dsolve.spsolve(Kmat, fext, use_umfpack=True)
	freac = Kreac.dot(d)
	return [d,freac]

def fesolve(Kmat,fext,ispc,vspc=None):
	n=np.shape(Kmat)[0]
	tr=0.0
	for i in xrange(n):
		tr += Kmat[i,i]
	beta=tr/(1.0*n)
	penfac=10000.0
	nspc=len(ispc)
	
	Kreac=Kmat[ispc,:].tocsr()
	if vspc==None:
		vspc=np.zeros((nspc,1))
	else: # inhomogenious bc
		for i in xrange(nspc):
			ii = ispc[i]
			iv = vspc[i]
			for jj in xrange(Kmat.shape[0]):
			     fext[jj] = fext[jj] - Kmat[jj,ii]*iv
			
	for i in xrange(nspc):
  	       ii = ispc[i]
	       Kmat[ii,ii] = penfac*beta
	       fext[ii] = penfac*beta*vspc[i]

	Kmat=Kmat.tocsr()
	d = dsolve.spsolve(Kmat, fext, use_umfpack=True)
	freac = Kreac.dot(d) 
	return [d,freac]
	
class FeSolver(object):
    
    def __init__(self,Kmat=None):
        self.Kmat = Kmat
        self.penfact = 10000.0
        
#    def Solve(self,Kmat,fext,spcs):
#        self.Kmat = Kmat
#        return self.Solve(fext,spcs)
        
    def Solve(self,fext,spcs):
        
        if ( self.Kmat==None ):
            print 'Kmat not set  in FeSolver'
            return
            
        n=np.shape(self.Kmat)[0]
        tr=0.0
        for i in xrange(n):
	    tr += self.Kmat[i,i]
	beta=tr/(1.0*n)
	#nspc=len(spcs)
	
	Kreac=self.Kmat.diagonal() 
	kappa = self.penfact*beta
	for ispc in spcs:
	   ii = ispc.gid
	   iv = ispc.dval
	   for jj in xrange(n):
	       fext[jj] = fext[jj] - self.Kmat[jj,ii]*iv
			
	for ispc in spcs:
  	   ii = ispc.gid
	   self.Kmat[ii,ii] = kappa
	   fext[ii] = kappa*ispc.dval

	self.Kmat=self.Kmat.tocsr()
	d = dsolve.spsolve(self.Kmat, fext, use_umfpack=True)
	
	for ii in xrange(n):
	    self.Kmat[ii,ii] = Kreac[ii]
	freac = self.Kmat.dot(d) 
	
	return [d,freac]
		
#------------
#ne=8
#kedim=6
#conn = np.array([[0,1,2,3,6,7],[2,3,8,9,6,7],[2,3,4,5,8,9],[4,5,10,11,8,9],\
#        [6,7,8,9,12,13],[8,9,14,15,12,13],[8,9,10,11,14,15],\
#        [10,11,16,17,14,15]],INDX_TYPE)
#ke = np.array([[7.4176,   -3.5714,   -1.9231,    1.6484,   -5.4945,    1.9231],\
#   [-3.5714,    7.4176,    1.9231,   -5.4945,    1.6484,   -1.9231],\
#   [-1.9231,    1.9231,    1.9231,         0,         0,   -1.9231],\
#   [ 1.6484,   -5.4945,         0,    5.4945,   -1.6484,         0],\
#   [-5.4945,    1.6484,         0,   -1.6484,    5.4945,         0],\
#   [ 1.9231,   -1.9231,   -1.9231,         0,         0,    1.9231]],FLOAT_TYPE)
#
#ne = 5
#kedim = 2
#ke = np.array([[100.0,-100.0],[-100.0,100.0]],FLOAT_TYPE)
#conn = np.array([[0,1],[1,2],[2,3],[3,4],[4,5]],INDX_TYPE)
#kdat = DelayAssm(ne*kedim**2)
#for e in conn:
#    kdat.AddKe(ke,e)

ke=np.array([[1,1,1,1],[1,1,1,1]],float)
kdat = DelayedAssm(16)
kdat.AddLocalMatrix(ke,np.array([0,1],int),np.array([0,1,2,3],int))
kdat.AddLocalMatrix(ke,np.array([1,2],int),np.array([2,3,4,5],int))