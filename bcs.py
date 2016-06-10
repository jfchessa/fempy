import numpy as np
import basic
import element as elem
import dofmap as dmap

import pdb

def compute_fext_tract(node,faceElem,trac,dofmap=dmap.FixedDofMap()):
    """Computes the external force for a constant traction"""
    nn = len(node)
    dpn = dmap.NumNodeDof()
    ndof=nn*dpn
    sdim=len(trac)
    fext = np.zeros( ndof, basic.FLOAT_TYPE )
    etype=None
    for e in faceElem:
        ecoord = node[e.Connectivity(),:]
        sctr = dofmap.Sctr(e.Connectivity(),range(dpn))
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
    
        # scatter fe into f
        fext[sctr] = fext[sctr] + fe
        
    return fext
    
    
class Force(object):
    
    def __init__(self,nid,f):
        
        self.nid = int(nid)
        self.fvect = np.array(f,basic.FLOAT_TYPE)
        self.gid = 6*self.nid + self.ldof
    
    def __repr__(self):
        return "Force on node "+str(self.nid)+" of "+str(self.fvect)

def SetGids(bcs,ndpn=6):
    for f in bcs:
        f.gid = ndpn * f.nid + f.ldof   

class PLoad(Force):
    
    def __init__(self,face,trac):
        self.face = face
        self.trac = trac
        
    def __repr__(self):
        return "Pressure load"

class SPC(object):
    
    def __init__(self,nid,ldof,dval=0.0):
        
        self.nid = int(nid)
        self.ldof = int(ldof)
        self.dval = float(dval)
        self.gid = 6*self.nid + self.ldof
    
    def __repr__(self):
        return "SPC on node "+str(self.nid)+" of "+str(self.dval)+" on ldof "+str(self.ldof)
  
class EssentialBCs(object):
    
    def __init__(self):
        self.spcs=[]
                            
    def AddSpcs(self,nids,ldofs,dvals=None,dofmap=dmap.FixedDofMap()):
        
        if ( dvals==None ):
            dvals = np.zeros( len(ldofs), basic.FLOAT_TYPE )
        
        for n in nids:
            for l in xrange( len(ldofs) ):
                self.spcs.append( SPC(n,ldofs[l],dvals[l]) )
                
        for s in self.spcs:
            s.gid = dofmap.GID( s.nid, s.ldof )
                
    def SPCs(self):
        return self.spcs
        
    def SetValues(self,U):
        for s in self.spcs:
            U[s.gid] = s.dval
          
class NaturalBCs(EssentialBCs):
    
    def __init__(self):
        self.forces=[]
        self.ploads=[]
        
    def AddPointForces(self,nids,fvct):
        
        for n in nids:
            self.forces.append( Force(n,fvct) )
                 
    def AddTraction(self,elem,trac):
        
        for e in elem:
            self.ploads.append( PLoad(e,trac) )
                                                                                    
    def Fext(self,node,dofmap=dmap.FixedDofMap()):
        
        nn = len(node)
        ndof = dofmap.GID(nn,0)
        fext = np.zeros( ndof, basic.FLOAT_TYPE )
        
        for f in self.forces:
            ldof = np.array( range(len(f.fvct)), basic.INDX_TYPE )
            sctr = dofmap.Sctr([f.nid],ldof)
            fext[sctr] = fext[sctr] + f.fvct
                
        for f in self.ploads:
            
            dpn = len(f.trac)
            fe = np.zeros( f.face.NumNodes()*dpn, basic.FLOAT_TYPE )
            ldof = np.array( range(dpn), basic.INDX_TYPE )
            sctr = dofmap.Sctr(f.face.Connectivity(),ldof)
            
            # compute fe
            etype=None
            e = f.face
            ecoord = node[e.Connectivity(),:]
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(e.Order())
                fedim = e.NumNodes()*dpn
            
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                N = e.N(qpt)
                jac = e.Jacobian(ecoord,qpt)
                for s in xrange(dpn):
                    fe[s:fedim:dpn] = fe[s:fedim:dpn] + N*f.trac[s]*jac*qwt
                    
            fext[sctr] = fext[sctr] + fe
            
        return fext
          