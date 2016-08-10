import numpy as np
from scipy import sparse
from scipy.sparse.linalg import dsolve

import pdb

###############################################################################
#
#                            -  F E M P Y   -  
#
###############################################################################


# numerical type definitions
FLOAT_TYPE  =  np.float64   # floats
#INDX_TYPE   =  np.uint32    # unsigned int
INDX_TYPE   =  np.int32     # ints for index things
        
        
#def scatter_matrix(ke,Kmat,sctr1,sctr2=None):
#	if ( sctr2==None ):
#		sctr2=sctr1
#	for i in xrange(len(sctr1)):
#		ii=sctr1[i]
#		for j in xrange(len(sctr2)):
#			jj=sctr2[j]
#			Kmat[ii,jj] += ke[i,j]

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
#def fesolve2(Kmat,fext,ispc,vspc=None):
#	n=np.shape(Kmat)[0]
#	tr=0.0
#	for i in xrange(n):
#		tr += Kmat[i,i]
#	beta=tr/(1.0*n)
#	nspc=len(ispc)
#	
#	Kreac=Kmat[ispc,:].tocsr()
#	if vspc==None:
#		vspc=np.zeros((nspc,1))
#	else: # inhomogenious bc
#		for i in xrange(nspc):
#			ii = ispc[i]
#			iv = vspc[i]
#			for jj in xrange(Kmat.shape[0]):
#			     fext[jj] = fext[jj] - Kmat[jj,ii]*iv
#			
#	for i in xrange(nspc):
#  	       ii = ispc[i]
#  	       for jj in xrange(n):
#  	           if ( not(Kmat[ii,jj]==0) ):
#	               Kmat[ii,jj] = 0.0
#	           if ( not(Kmat[jj,ii]==0) ):
#	               Kmat[jj,ii] = 0.0
#	       Kmat[ii,ii] = beta
#	       fext[ii] = beta*vspc[i]
#
#	Kmat=Kmat.tocsr()
#	d = dsolve.spsolve(Kmat, fext, use_umfpack=True)
#	freac = Kreac.dot(d)
#	return [d,freac]

#def fesolve(Kmat,fext,ispc,vspc=None):
#	n=np.shape(Kmat)[0]
#	tr=0.0
#	for i in xrange(n):
#		tr += Kmat[i,i]
#	beta=tr/(1.0*n)
#	penfac=10000.0
#	nspc=len(ispc)
#	
#	Kreac=Kmat[ispc,:].tocsr()
#	if vspc==None:
#		vspc=np.zeros((nspc,1))
#	else: # inhomogenious bc
#		for i in xrange(nspc):
#			ii = ispc[i]
#			iv = vspc[i]
#			for jj in xrange(Kmat.shape[0]):
#			     fext[jj] = fext[jj] - Kmat[jj,ii]*iv
#			
#	for i in xrange(nspc):
#  	       ii = ispc[i]
#	       Kmat[ii,ii] = penfac*beta
#	       fext[ii] = penfac*beta*vspc[i]
#
#	Kmat=Kmat.tocsr()
#	d = dsolve.spsolve(Kmat, fext, use_umfpack=True)
#	freac = Kreac.dot(d) 
#	return [d,freac]
	
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
    
    def __init__(self,*args):
        self.nodes = {}

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
        
   
class DofMap(object):
    
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
                    
    def ActivateDofs(self,elements,ldofs):
        for eid, e in elements.iteritems():
            for n in e.conn:
                
                if not n in self.gid:
                    self.gid[n] = {}
                
                for s in ldofs:
                    if not s in self.gid[n]:
                        self.gid[n][s] = self.ndof
                        self.ndof += 1
                        
    def Renumber(self):
        I=0
        for i, gdofs in self.gid.iteritems():
            n = len(gdofs)
            self.gid[i] = range(I,I+n)
            I += n
              
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
        return np.array( self.gid[n].keys() )
            
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
            
    def __repr__(self):
        s='EssentialBCs\nnid: ldof values\n'
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
        for e in faceelem:
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
        
        for n, ldof in self.iteritems():
            for s, v in ldof.iteritems():
                rhs[ dofmap.GID(n,s) ] += v
        
# ---------------------------------------------------------------------------
#
#              Quadrature Rules
#
# ----------------------------------------------------------------------------

def quadrature_gauss1d(numpt):
    
    if ( numpt==1 ):
        quadpoint = np.array( [0.000000000000000], dtype =  FLOAT_TYPE )
        quadweight = np.array( [2.000000000000000], dtype =  FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 2 ):
        quadpoint = np.array( [0.577350269189626, -0.577350269189626], dtype =  FLOAT_TYPE )
        quadweight = np.array( [1.0, 1.0 ], dtype =  FLOAT_TYPE )  
        return [quadpoint,quadweight]
        
    if ( numpt == 3 ):
        quadpoint = np.array( [0.774596669241483, -0.774596669241483, 0.000000000000000], dtype =  FLOAT_TYPE )
        
        quadweight = np.array( [0.555555555555556, 0.555555555555556, 0.888888888888889], dtype =  FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 4 ):
        quadpoint = np.array( [0.861134311594053,-0.861134311594053, \
                0.339981043584856,-0.339981043584856], dtype =  FLOAT_TYPE )
        
        quadweight = np.array( [0.347854845137454, 0.347854845137454,\
                0.652145154862546, 0.652145154862546], dtype =  FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 5 ):
        quadpoint = np.array( [0.906179845938664,-0.906179845938664,0.538469310105683,\
             -0.538469310105683, 0.000000000000000], dtype =  FLOAT_TYPE )
        
        quadweight = np.array( [0.236926885056189, 0.236926885056189, 0.478628670499366,\
                0.478628670499366, 0.568888888888889], dtype =  FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 6 ):
        quadpoint = np.array( [0.932469514203152,-0.932469514203152,\
            0.661209386466265,-0.661209386466265,0.238619186003152,\
            -0.238619186003152], dtype =  FLOAT_TYPE )
        
        quadweight = np.array( [0.171324492379170,0.171324492379170,0.360761573048139,\
                0.360761573048139,0.467913934572691, 0.467913934572691], dtype =  FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 7 ):
        quadpoint = np.array( [0.949107912342759,-0.949107912342759,0.741531185599394,\
                            -0.741531185599394,0.405845151377397,\
                            -0.405845151377397,0.000000000000000], dtype =  FLOAT_TYPE )
    
        quadweight = np.array( [0.129484966168870,0.129484966168870,0.279705391489277,\
                0.279705391489277,0.381830050505119,0.381830050505119,\
                0.417959183673469 ], dtype =  FLOAT_TYPE )  
        return [quadpoint,quadweight]
        
    if ( numpt == 8 ):
        quadpoint = np.array( [0.960289856497536,-0.960289856497536,0.796666477413627,\
         -0.796666477413627,0.525532409916329,-0.525532409916329,0.183434642495650,\
         -0.183434642495650], dtype =  FLOAT_TYPE )
        
        quadweight = np.array( [0.101228536290376, 0.101228536290376, 0.222381034453374,\
            0.222381034453374, 0.313706645877887, 0.313706645877887,\
            0.362683783378362, 0.362683783378362 ], dtype =  FLOAT_TYPE )  
        return [quadpoint,quadweight]

#def quadrature_gauss( sdim, p ):
#    
#    nn = int(np.ceil(0.5*(p+1)))
#    numpt = nn ** sdim
#    if ( sdim == 1 ):
#        return quadrature_gauss1d( numpt )
#
#    qrule = quadrature_gauss1d(nn)
#
#    quadweight = np.zeros( numpt, dtype =  FLOAT_TYPE )
#    quadpoint = np.zeros( (numpt,sdim), dtype =  FLOAT_TYPE ) 
#    
#    for s in xrange(sdim):
#        quadweight[s*nn:(s+1)*nn]  = qrule[1]
#        quadpoint[s*nn:(s+1)*nn,0] = qrule[0]
#        quadpoint[s:numpt:nn,1] = qrule[0]
#        
#    return [quadpoint,quadweight]
            
def compound_quadrature(p1,w1,p2,w2,p3=None,w3=None):
    
    if ( p3==None ): 
        
        npts = w1.size * w2.size
        sdim = 2 
        qwts = np.zeros( npts, dtype =  FLOAT_TYPE )
        qpts = np.zeros((npts,sdim), dtype =  FLOAT_TYPE )
        
        n=0
        for j in xrange(w2.size):        
            for i in xrange(w1.size):  
                qpts[n,:]  =  [ p1[i], p2[j] ]  
                qwts[n] = w1[i]*w2[j]
                n = n+1
        
    else:
        
        npts = w1.size * w2.size * w3.size
        sdim = 3
        qwts = np.zeros( npts, dtype =  FLOAT_TYPE )
        qpts = np.zeros((npts,sdim), dtype =  FLOAT_TYPE )
        
        n=0
        for k in xrange(w3.size):
            for j in xrange(w2.size):        
                for i in xrange(w1.size):  
                    qpts[n,:]  =  [ p1[i], p2[j], p3[k] ]  
                    qwts[n] = w1[i]*w2[j]*w3[k]
                    n = n+1
                          
    return [qpts,qwts]
                                                                      
def quadrature_simplex(sdim,quadorder):
    sixth = 1.0/6
    
    if ( sdim == 3 ):  # tetrahedra
#      if ( quadorder ~= 1 &  quadorder ~= 2 &  quadorder ~= 3  ) 
#        % check for valid quadrature order
#        disp('Incorect quadrature order for QUADRATURE_SIMPLEX');
#        quadorder = 1;
#       end
        
        if  ( quadorder == 1 ):
            quadpoint = np.array([[ 0.25, 0.25, 0.25 ]],dtype= FLOAT_TYPE)
            quadweight = np.array([ 1.0 ],dtype= FLOAT_TYPE)
            return [quadpoint,sixth*quadweight]
        
        if ( quadorder == 2 ): 
            quadpoint = np.array([[ 0.58541020,  0.13819660,  0.13819660],\
                      [ 0.13819660,  0.58541020,  0.13819660],\
                      [ 0.13819660,  0.13819660,  0.58541020],\
                      [ 0.13819660,  0.13819660,  0.13819660]],dtype= FLOAT_TYPE)
            quadweight = np.array([.25,.25,.25,.25],dtype= FLOAT_TYPE)
            return [quadpoint,sixth*quadweight]
        
        if ( quadorder == 3 ):
            quadpoint = np.array( [[ 0.25,  0.25,  0.25],\
                     [ 0.5,   sixth,   sixth ],\
                     [ sixth,   0.5,   sixth ],\
                     [ sixth,   sixth,   0.5 ],\
                     [ sixth,   sixth,   sixth ]], dtype= FLOAT_TYPE )
            quadweight = np.array( [-0.8, .45, .45, .45, .45], dtype= FLOAT_TYPE )
            return [quadpoint,sixth*quadweight]
        
         
    else:  # TRIANGLES
      
#      if ( quadorder > 7 ) % check for valid quadrature order
#        disp('Quadrature order too high for QUADRATURE_SIMPLEX');
#        quadorder = 1;
#      end
      
        if ( quadorder <= 1 ):  # set quad points and quadweights
            quadpoint = np.array([[ 0.3333333333333, 0.3333333333333 ]], dtype= FLOAT_TYPE )
            quadweight =np.array([ 0.5 ], dtype= FLOAT_TYPE )
            return [quadpoint,quadweight]  
        
        if ( quadorder == 2 ): 
            quadweight = sixth*np.ones( 3, dtype= FLOAT_TYPE )
            quadpoint = np.array( [[ 0.1666666666667, 0.1666666666667 ],\
                            [ 0.6666666666667, 0.1666666666667 ],\
                            [ 0.1666666666667, 0.6666666666667 ]], dtype= FLOAT_TYPE ) 
            return [quadpoint,quadweight]  
        
        if ( quadorder <= 5 ): 
        
            quadpoint = np.array( [[ 0.1012865073235, 0.1012865073235 ],\
                        [ 0.7974269853531, 0.1012865073235 ],\
                        [ 0.1012865073235, 0.7974269853531 ],\
                        [ 0.4701420641051, 0.0597158717898 ],\
                        [ 0.4701420641051, 0.4701420641051 ],\
                        [ 0.0597158717898, 0.4701420641051 ],\
                        [ 0.3333333333333, 0.3333333333333 ]], dtype= FLOAT_TYPE )
        
            quadweight = np.array( [0.1259391805448, 0.1259391805448, 0.1259391805448,\
                    0.1323941527888, 0.1323941527885,  0.1323941527885, \
                    0.2250000000000 ], dtype= FLOAT_TYPE )
                    
            return [quadpoint,0.5*quadweight]    
        
        quadpoint = np.array( [[ 0.0651301029022, 0.0651301029022 ],\
                    [ 0.8697397941956, 0.0651301029022 ],\
                    [ 0.0651301029022, 0.8697397941956 ],\
                    [ 0.3128654960049, 0.0486903154253 ],\
                    [ 0.6384441885698, 0.3128654960049 ],\
                    [ 0.0486903154253, 0.6384441885698 ],\
                    [ 0.6384441885698, 0.0486903154253 ],\
                    [ 0.3128654960049, 0.6384441885698 ],\
                    [ 0.0486903154253, 0.3128654960049 ],\
                    [ 0.2603459660790, 0.2603459660790 ],\
                    [ 0.4793080678419, 0.2603459660790 ],\
                    [ 0.2603459660790, 0.4793080678419 ],\
                    [ 0.3333333333333, 0.3333333333333 ]], dtype= FLOAT_TYPE )
        
        quadweight = np.array( [ 0.0533472356088, 0.0533472356088, 0.0533472356088,\
                    0.0771137608903, 0.0771137608903, 0.0771137608903, 0.0771137608903,\
                    0.0771137608903, 0.0771137608903, 0.1756152576332, 0.1756152576332, \
                    0.1756152576332, -0.1495700444677 ], dtype= FLOAT_TYPE )
                     
        return [quadpoint,0.5*quadweight]  
        
    
                
#  ---------------------------------------------------------------------------
#
#             Element classes and support functions
#
#  ---------------------------------------------------------------------------

def sctr_array(econn,ndofn):
    nn = len(econn)
    ss = nn*ndofn
    sctr = np.ndarray( ss,  INDX_TYPE )
    s0 = ndofn*np.array(econn, INDX_TYPE)
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
    
    jmat = np.zeros( (sdim,sdim), dtype= FLOAT_TYPE )
    
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
            jmat[:,1] = np.cross( np.array([0,0,1],dtype= FLOAT_TYPE), jmat[:,0] )
            if ( np.linalg.norm(jmat[:,1])==0 ):
                jmat[:,1] = np.cross( np.array([0,1,0],dtype= FLOAT_TYPE), jmat[:,0] )
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
    bmat=np.zeros([3,2*nn], dtype= FLOAT_TYPE)
    bmat[0,0:2*nn:2] = dndx.transpose()[0,:]
    bmat[1,1:2*nn:2] = dndx.transpose()[1,:]
    bmat[2,0:2*nn:2] = dndx.transpose()[1,:]
    bmat[2,1:2*nn:2] = dndx.transpose()[0,:]
    return bmat
    
    
def form_bmat_3d(dndx):
    nn=dndx.shape[0]
    bmat=np.zeros([6,3*nn], dtype= FLOAT_TYPE)
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
        n = int(math.ceil(0.5*(p+1)))
        qpt, qwt =  quadrature_gauss1d(n)
        return [ np.reshape(qpt,(len(qpt),1)), qwt ]
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        basis = (0.5)*np.array( [(1-pt[0]), (1+pt[0]) ], dtype= FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        dbasis = np.array([[ -0.5 ],[ 0.5 ]], dtype= FLOAT_TYPE )
        return dbasis

    def dNdx(self,coords,pt=None,sdim=None):
        
        if ( coords.ndim == 1 ): 
            coords = np.reshape( coords, (coords.size,1) )
        
        if ( sdim == None ):
            sdim = coords.shape[1]
        
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
              
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
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
            
        dNdx = self.dNdx(coords,pt,sdim)
        
        if ( self.ElemDim() == sdim ):
            return dNdx
            
        bmat = np.zeros( (1,sdim*self.NumNodes()), dtype= FLOAT_TYPE )
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
#  Quadrature Rules
#
#


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
        return  quadrature_simplex(self.ElemDim(),p)
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        basis = np.array( [1-pt[0]-pt[1], pt[0], pt[1] ], dtype= FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        dbasis = np.array([ [-1,-1], [1,0], [0,1] ], dtype= FLOAT_TYPE )
        return dbasis
        
    def BMat(self,coords,pt=None,sdim=None):
        
        if ( coords.ndim == 1 ): 
            coords = np.reshape( coords, (coords.size,1) )
        
        if ( sdim == None ):
            sdim = coords.shape[1]
        
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
            
        ipm = self.dNdx(coords,pt,sdim)
        
        if ( self.ElemDim() == sdim ):
            return [ form_bmat_2d(ipm[0]), ipm[1] ]     
 
        bmat = np.zeros( (6,sdim*self.NumNodes()), dtype= FLOAT_TYPE )
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
        n = int(math.ceil(0.5*(p+1)))
        qr =  quadrature_gauss1d(n)
        return  compound_quadrature(qr[0],qr[1],qr[0],qr[1])
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        basis = (0.25)*np.array( [(1-pt[0])*(1-pt[1]), \
                     (1+pt[0])*(1-pt[1]), \
                     (1+pt[0])*(1+pt[1]), \
                     (1-pt[0])*(1+pt[1]) ], dtype= FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        dbasis = (0.25)*np.array([ [ -(1-pt[1]), -(1-pt[0]) ],\
                                   [  (1-pt[1]), -(1+pt[0]) ],\
                                   [  (1+pt[1]),  (1+pt[0]) ],\
                                   [ -(1+pt[1]),  (1-pt[0]) ]], dtype= FLOAT_TYPE )
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
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        basis = np.array( [1-pt[0]-pt[1]-pt[2], pt[0], pt[1], pt[2] ], dtype= FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        dbasis = np.array([ [-1,-1,-1], [1,0,0], [0,1,0], [0,0,1] ], dtype= FLOAT_TYPE )
        return dbasis
        
    def BMat(self,coords,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
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
        n = int(math.ceil(0.5*(p+1)))
        qr =  quadrature_gauss1d(n)
        return  compound_quadrature(qr[0],qr[1],qr[0],qr[1],qr[0],qr[1])
        
    def N(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        basis = .125*np.array( [(1-pt[0])*(1-pt[1])*(1-pt[2]), \
                     (1+pt[0])*(1-pt[1])*(1-pt[2]), \
                     (1+pt[0])*(1+pt[1])*(1-pt[2]), \
                     (1-pt[0])*(1+pt[1])*(1-pt[2]), \
                     (1-pt[0])*(1-pt[1])*(1+pt[2]), \
                     (1+pt[0])*(1-pt[1])*(1+pt[2]), \
                     (1+pt[0])*(1+pt[1])*(1+pt[2]), \
                     (1-pt[0])*(1+pt[1])*(1+pt[2]) ], dtype= FLOAT_TYPE )
        return basis
        
    def dNdxi(self,pt=None):
        if ( pt==None ):
            pt = np.zeros( self.ElemDim(), dtype= FLOAT_TYPE )
        dbasis = 0.125*np.array([ [-(1-pt[1])*(1-pt[2]), -(1-pt[0])*(1-pt[2]), -(1-pt[0])*(1-pt[1]) ],\
                                   [  (1-pt[1])*(1-pt[2]), -(1+pt[0])*(1-pt[2]), -(1+pt[0])*(1-pt[1]) ],\
                                   [  (1+pt[1])*(1-pt[2]),  (1+pt[0])*(1-pt[2]), -(1+pt[0])*(1+pt[1]) ],\
                                   [ -(1+pt[1])*(1-pt[2]),  (1-pt[0])*(1-pt[2]), -(1-pt[0])*(1+pt[1]) ],\
                                   [ -(1-pt[1])*(1+pt[2]), -(1-pt[0])*(1+pt[2]),  (1-pt[0])*(1-pt[1]) ],\
                                   [  (1-pt[1])*(1+pt[2]), -(1+pt[0])*(1+pt[2]),  (1+pt[0])*(1-pt[1]) ],\
                                   [  (1+pt[1])*(1+pt[2]),  (1+pt[0])*(1+pt[2]),  (1+pt[0])*(1+pt[1]) ],\
                                   [ -(1+pt[1])*(1+pt[2]),  (1-pt[0])*(1+pt[2]),  (1-pt[0])*(1+pt[1]) ] ], dtype= FLOAT_TYPE )
        return dbasis
 
    	
class ElementArray(dict):
    
    def __init__(self,*args):
        self.elems = {}

 
    def __repr__(self):
        s='Element Array\n'
        for eid, e in self.iteritems():
            s=s+str(eid)+': ' +str(e)+'\n'
        return s
               	    	    	                         	   	    	     
               	   	    	       	    	     
               	   	    	       	    	       	    	       	    	     
               	    	  
###########################################################################
pwidth=10
pheight=20

node = NodeArray()
n=0
nnx=3
nny=5
for y in np.linspace(0,pheight,nny):
    for x in np.linspace(0,pwidth,nnx):
        node.AddNode(n,[x,y,0.0])
        n += 1

element = ElementArray()
prop=0
conn=np.array([0, 1, 6, 5])
e=0
for j in xrange(nny-1):
    for i in xrange(nnx-1): 
        element[e] = ElemQuad4(conn,prop) 
        conn += 1
        e += 1
    conn += 1
    
dofmap = DofMap()
dofmap.ActivateDofs(element,[0,1])