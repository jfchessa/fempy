import numpy as np

#import scipy as sp
from scipy import sparse
from scipy.sparse.linalg import dsolve

#import pdb

# numerical type definitions
FLOAT_TYPE  =  np.float64   # floats
#INDX_TYPE   =  np.uint32    # unsigned int
INDX_TYPE   =  np.int32     # ints for index things
        
        
def scatter_matrix(ke,Kmat,sctr1,sctr2=None):
    """scatter_matrix(ke,Kmat,sctr1,sctr2=None)"""
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
        for eid, e in element.iteritems():
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

		
#class NodeArray(object):
#    
#    def __init__(self,cap=None):
#        self.sinc =100
#        if cap==None:
#            cap=self.sinc
#        
#        self.ncoord = np.zeros((cap,3),dtype=np.float)
#        self.nid = {}
#        
#        self.nn = 0
#         
#    def __repr__(self):
#        s='Node Array\n'
#        for k, v in self.nid.iteritems():
#            s=s+str(k)+': ' +str(self.ncoord[v])+'\n'
#        return s
#                
#    def Capacity(self):
#        return len(self.ncoord)   
#        
#    def IncreaseCapacity(self,cap):
#        if cap<=self.Capacity():
#            return
#        
#        tmp = self.ncoord[range(self.nn)]
#        self.ncoord = np.zeros((cap,3),dtype=np.float)
#        self.ncoord[range(self.nn)] = tmp     
#            
#    def AddNode(self,n,coord):
#        
#        if self.nn >= self.Capacity():
#            self.IncreaseCapacity(self.Capacity()+self.sinc)
#            
#        self.nid[n] = self.nn
#        self.ncoord[self.nn] = coord       
#            
#        self.nn = self.nn + 1
#        
#    def __getitem__(self,i):
#        try:
#            return self.CoordMat(i)
#        except:
#            return self.CoordMat([i])
#                
#    def CoordMat(self,conn):
#        Is = np.zeros(len(conn),dtype=int)
#        for i in xrange(len(conn)):
#            try:
#                I = self.nid[conn[i]]
#            except KeyError:
#                print 'invalid node id'
#                return 
#            Is[i]= I
#       
#        return self.ncoord[Is]	
  

class Point(np.ndarray):
 
    def __new__( subtype, coord=[] ): 
        obj = np.ndarray.__new__( subtype, 3, float )
        return obj
        
    def __init__( self, coord=[] ):
        i=0
        for c in coord:
            self[i] = c
            i += 1
            if i > 2:
                return
        while i<3:
            self[i] = 0.0
            i += 1
            
    def __repr__( self ):
        return '['+str(self[0])+' '+str(self[1])+' '+str(self[2])+']'
        
    def x(self):
        return self[0]
        
    def y(self):
        return self[1]
        
    def z(self):
        return self[2]
        
        
class NodeArray(dict):
    
    def __init__(self,**kwargs):
        self.nodes = {}
        if 'narray' in kwargs:
            self.AddNodeArray( kwargs['narray'] )

    def __repr__(self):
        s='Node Array\n'
        for nid, node in self.iteritems():
            s=s+str(nid)+': ' +str(node)+'\n'
        return s
    
    def CoordMat(self,nids):
        cm = np.ndarray((len(nids),3),dtype=float)
        ii=0
        for i in nids:
            cm[ii]= dict.__getitem__(self,i)
            ii += 1
        return cm           	    	    	     
     
    def __getitem__(self,i):
        try:
            return self.CoordMat(i)
        except:
            return self.CoordMat([i])
                          	    	    	
    def __setitem__(self,i,node):
        dict.__setitem__(self,i,Point(node))
                
    def AddNode(self,n,node):
        self.__setitem__(n,Point(node))   
        
    def AddNodeArray(self,ncoord,nids=None):
        n = len(ncoord)
        if nids==None:
            nids=range(n)
            
        for i in range(n):
            self.AddNode( nids[i], ncoord[i] )   
        
    def NumNodes(self):
        return len(self)  
        
    def SDIM(self):
        return 3
        
    def ContinousNIDMap(self):
        """returns a map of node id to the location in the map ordering. 
        This is useful if you need to address the nodes in zero offset and
        continous numberng form."""
        nmap = dict()
        
        n=0
        for nid in self.iteritems():
            nmap[nid] = n
            n += 1
            
        return nmap
        
        
    def XYZRange(self):
        """returns the  x, y and z range of the nodes in the array"""
        
        for nid, node in self.iteritems():
            xmin=node[0]
            xmax=node[0]
            ymin=node[1]
            ymax=node[1]
            zmin=node[2]
            zmax=node[2]
            break
        
        for nid, node in self.iteritems():
            if xmin>node[0]:
                xmin=node[0]
            elif xmax<node[0]:
                xmax=node[0]
                
            if ymin>node[1]:
                ymin=node[1]
            elif ymax<node[1]:
                ymax=node[1]
                
            if zmin>node[2]:
                zmin=node[2]
            elif zmax<node[2]:
                zmax=node[2]
        
        return [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
        
    def Corners(self):
        """returns the corner verticies of the node array"""
        r = self.XYZRange()
        corners = np.array( [[ r[0][0], r[1][0], r[2][0] ],
                             [ r[0][1], r[1][0], r[2][0] ],
                             [ r[0][1], r[1][1], r[2][0] ],
                             [ r[0][0], r[1][1], r[2][0] ],
                             [ r[0][0], r[1][0], r[2][1] ],
                             [ r[0][1], r[1][0], r[2][1] ],
                             [ r[0][1], r[1][1], r[2][1] ],
                             [ r[0][0], r[1][1], r[2][1] ] ] )
        return corners
   
class DofMapFixed(object):
    
    def __init__(self,ndpn=6):
        self.ndpn=ndpn
 
    def __repr__(self):
        s='Dof Map Fixed, num dof per node: '+str(self.ndpn)+' \n'
        return s
                
    def GID(self,nid,lid=0):
        return self.ndpn * nid + lid
        
    def LDOFs(self,nid=0):    
        return np.array(range(self.ndpn),dtype=int)
            
    def Sctr(self,conn,ldofs=None):
        if ( ldofs==None ):
            ldofs = self.LDOFs()
            
        sctr = np.ndarray( len(conn)*len(ldofs), INDX_TYPE )
        ii = 0
        for nid in conn:
            for s in ldofs:
                sctr[ii] = self.GID(nid,s)
                ii = ii + 1
        return sctr   
        
class DofMap(DofMapFixed):
    
    def __init__(self,ndpn={}):
        self.gid = {}
        self.ndof=0
        for i, n in ndpn.iteritems():
            self.gid[i] = range(self.ndof,self.ndof+n)
            self.ndof += n
 
    def __repr__(self):
        s='Dof Map\n'
        for n, gdofs in self.gid.iteritems():
            s = s +str(n)+": "+str(gdofs)+'\n'
        return s
       
    #def __getitem__(self,n):
    #    try:
    #        s=[]
    #        for i in n:
    #            try:
    #                s=s.append(self.gid[i])
    #            except KeyError:
    #                s=s.append({})
    #        return s
    #    except TypeError:
    #        try:
    #            return self.gid[n]
    #        except KeyError:
    #            return {}
         
    #def ActivateDofs(self,elements,ndof):
    #    for eid, e in elements.iteritems():
    #        for n in e.conn:
    #            
    #            if n in self.gid:               # node already has active gdofs
    #                nn = ndof-self.NumNodeDof(n)
    #                if nn > 0:                  # but we need to add more
    #                    self.gid[n].append( range(self.ndof,self.ndof+nn) )
    #                    self.ndoe += nn
    #            else:                          # node not active 
    #                self.gid[n] = range(self.ndof,self.ndof+ndof)
    #                self.ndof += ndof
                    
    def ActivateDofs(self,elements,ldofs,renum=True):
        """  ActivateDofs(self,elements,ldofs,renum=True) """
        for eid, e in elements.iteritems():
            for n in e.conn:
                
                if not n in self.gid:
                    self.gid[n] = {}
                
                for s in ldofs:
                    if not s in self.gid[n]:
                        self.gid[n][s] = self.ndof
                        self.ndof += 1
        if renum:
            self.Renumber()
                        
    def Renumber(self):
        I=0
        for i, gdofs in self.gid.iteritems():
            for s, v in gdofs.iteritems():
                self.gid[i][s] = I
                I += 1
              
    def GID(self,nid,lid=0):
        try:
            return self.gid[nid][lid]
        except KeyError:
            print 'error in DofMap.GID, lid '+str(lid)+' not defined for node '+str(nid)
            return KeyError
        
    def NumNodeDof(self,nid=0):
        """Number of dof per node"""
        return len(self.gid[nid])
        
    def NumDof(self):
        """Total number of dofs"""
        return self.ndof
        
    def LDOFs(self,nid=0):    
        return np.array( self.gid[nid].keys() )
            
    def Sctr(self,conn,ldofs=None):
        if ( ldofs==None ):
            ldofs = self.LDOFs()
            
        sctr = np.ndarray( len(conn)*len(ldofs), INDX_TYPE )
        ii = 0
        for nid in conn:
            for s in ldofs:
                sctr[ii] = self.GID(nid,s)
                ii = ii + 1
        return sctr


class EssentialBCs(dict):
    """Class to define single point constraints and other constraints"""        
    def __repr__(self):
        s='EssentialBCs\nnid: ldof values\n'
        for n, ldof in self.iteritems():
            s=s+str(n)+': '+str(ldof)+'\n'
        return s
                
    def AddPointValues(self,nids,ldofs,vals=None):
        """AddPointValues(self,nids,ldofs,vals=None)"""
        try:
            nspc = len(ldofs)
        except TypeError:
            nspc = 1
            ldofs=[ldofs]
        
        if vals==None:
            vals=np.zeros(nspc)
    
        if not len(vals)==nspc:
            vals = np.ones(nspc)*vals[0]
            
        for n in nids:
            if not n in self:
                dict.__setitem__(self,n,{})
            for s in xrange(nspc):
                (dict.__getitem__(self,n))[ ldofs[s] ] = vals[s]    
                
    def GetIFIX(self,dofmap):
        
        ifix = []
        ival = []
        for n, ldof in self.iteritems():
            for s, v in ldof.iteritems():
                ifix.append( dofmap.GID(n,s) )
                ival.append( v )
        
        return [ ifix, ival ]
 
class NaturalBCs(dict):
            
    def __repr__(self):
        s='NaturalBCs\nnid: ldof values\n'
        for n, ldof in self.iteritems():
            s=s+str(n)+': '+str(ldof)+'\n'
        return s
                
    def AddPointValues(self,nids,ldofs,vals=None):
        
        nspc = len(ldofs)
        
        if vals==None:
            vals=np.zeros(nspc)
    
        if not len(vals)==nspc:
            vals = np.ones(nspc)*vals[0]
            
        for n in nids:
            if not n in self:
                dict.__setitem__(self,n,{})
            for s in xrange(nspc):
                if ldofs[s] in self[n]: 
                    (dict.__getitem__(self,n))[ ldofs[s] ] += vals[s] 
                else:  
                    (dict.__getitem__(self,n))[ ldofs[s] ] = vals[s]  
              
    def AddFaceTraction(self,node,faceelem,itrac,trac):
        dpn = len(trac)
        sdim=len(trac)
        etype=None
        
        try:
            eiter = faceelem.values()
        except AttributeError:
            eiter = faceelem
        
        for e in eiter:
            ecoord = node[e.Connectivity()]
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                fedim = e.NumNodes()*dpn
                
            fe = np.zeros(fedim)
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                N = e.N(qpt)
                jac = e.Jacobian(ecoord,qpt)
                for s in xrange(sdim):
                    fe[s:fedim:sdim] = fe[s:fedim:sdim] + N*trac[s]*jac*qwt
        
            ii=0
            for i in e.Connectivity():
                self.AddPointValues([i],itrac,fe[ii:ii+sdim])
                ii += sdim
          
    def AddRHS(self,dofmap,rhs):
        """NatrualBCs.AddRHS(dofmap,rhs) : adds natural bcs valeus to rhs vector"""
        for n, ldof in self.iteritems():
            for s, v in ldof.iteritems():
                rhs[ dofmap.GID(n,s) ] += v
        
                
