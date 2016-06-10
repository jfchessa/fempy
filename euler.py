import numpy as np
import basic
import element as elem  
import dofmap as dmap

import pdb
              
# -------------------------------------------------------------------------   
class PolytropicGas(object):
    
    def __init__(self,cp=1.005,cv=.718):
        self.gamma = cp/cv
        self.cv = cv
        
    def GetP(self,rho,v,e):
        ke=0.0
        for vi in v:
            ke = ke + vi*vi
        return (self.gamma-1)*rho*(e-0.5*ke)
        
    def GetT(self,rho,v,e):
        ke=0.0
        for vi in v:
            ke = ke + vi*vi
        return (self.cv)*(e-0.5*ke)
        
    def GetS(self,rho,v,e):
        sdim = len(v)
        return np.zeros( (sdim,sdim), basic.FLOAT_TYPE )

class EulerDofMap(dmap.FixedDofMap):
    
    def __init__(self,nn,sdim=3):
        self.dpn = sdim + 2
        self.nn = nn
      
    def GID(self,nid,lid):
        if ( lid==0 ):
            return nid
        elif ( lid<self.dpn ):
            return self.sdim*nid + lid + self.nn
        else:
            return (1+self.sdim)*self.nn + nid     
    
class CompressibleDofMap(dmap.FixedDofMap):
    
    def __init__(self,nn,sdim=3):
        self.dpn = sdim + 1
        self.nn = nn
      
    def GID(self,nid,lid):
        if ( lid==0 ):
            return nid
        else:
            return self.sdim*nid + lid + self.nn
    
class CBSMethod(object):
    
    def __init( self, fluid=PolytropicGas(), sdim=3 ):
        self.sdim = sdim
        self.fluid = fluid
        self.dpn = self.sdim+2
        
        self.DofMap = None
        self.U  = None
    
    def Initialize(self,node,element):
        nn = len(node)
        ndof = self.dpn*nn
        
        self.DofMap = EulerDofMap(nn,self.sdim)
        self.Fa = np.zeros(ndof,basic.FLOAT_TYPE)
        self.Fv = np.zeros(ndof,basic.FLOAT_TYPE)
        self.M  = np.zeros(ndof,basic.FLOAT_TYPE)
        self.U  = np.zeros(ndof,basic.FLOAT_TYPE)
    
    def Step1(self,node,element):
        """Compute the advection term"""
        
        if ( self.Fa == None ):
            self.Initialize(node,element)
                
        etype = None
        for e in element:
            
            if ( not(etype==e.Type()) ):
                nn = e.NumNodes()
                
            fe = np.zeros(self.dpn*nn,basic.FLOAT_TYPE)
            sctr = self.DofMap.Sctr(e.conn)
            
            