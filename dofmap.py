import numpy as np
import basic

    
class DofMap(object):
    
    def __init__(self,ndpn):
        nn = len(ndpn)
        self.gidn=np.array(nn+1,basic.INDX_TYPE)
        self.gidn[0]=0
        for i in xrange(nn):
            self.gidn[i+1] = self.gidn[i] + ndpn[i]
            
    def GID(self,nid,lid=0):
        return self.ndpn[nid] + lid
        
    def NumNodeDof(self,nid=0):
        """Number of dof per node"""
        return self.ndpn[nid]
        
    def LDOFs(self,nid=0):    
        return np.array( range(self.NumNodeDof()), int )
            
    def Sctr(self,conn,ldofs=None):
        if ( ldofs==None ):
            ldofs = self.LDOFs()
            
        sctr = np.ndarray( len(conn)*len(ldofs), basic.INDX_TYPE )
        ii = 0
        for nid in conn:
            for s in ldofs:
                sctr[ii] = self.GID(nid,s)
                ii = ii + 1
        return sctr

class FixedDofMap(DofMap):
    
    def __init__(self,ndpn=6):
        self.ndpn = ndpn
        
    def GID(self,nid,lid):
        return self.ndpn * nid + lid
        
    def NumNodeDof(self,nid=0):
        return self.ndpn
        
class VariDofMap(DofMap):
    
    def __init__(self,nids):
        
        self.dmap = { }
        self.fixed = False
        self.ndof=0
        
        if type(nids)==dict:
            for k, n in nids.iteritems():
                self.dmap[k] = []
        else:            
            for n in nids:
                self.dmap[n] = []
            
    def ActivateDofs(self,elements,ndof):
        for e in elements:
            for n in e.conn:
                self.dmap[n] = [-1]*ndof
                self.ndof+=ndof
                
    def Fixed(self):
        return self.fixed
        
    def NumDof(self):
        return self.ndof
                
    def NumNodeDof(self,nid):
        """Number of dof per node"""
        return len(self.dmap[nid])
    
    def LDOFs(self,nid):    
        return self.dmap[nid]
            
    def FixDofs(self):
        if self.fixed:
            return
        n=0
        for key, val in self.dmap.iteritems():
            ndof  = len(val)
            self.dmap[key] = range(n,n+ndof)
            n += ndof
        self.fixed = True      
        
    def GID(self,nid,lid):
        if self.fixed:
            return self.dmap[nid][lid]
        return -1
