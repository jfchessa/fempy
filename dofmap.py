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
        
    def LDOFs(self):    
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
        

        