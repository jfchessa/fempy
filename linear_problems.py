import numpy as np
import scipy as sp
from scipy import sparse

import basic
import material as matl
import element as elem  
import dofmap as dmap

import pdb

        
                                    
# -------------------------------------------------------------------------        
class ProblemPlaneStress(object):
    
    def __init__(self):
        self.sdim = 2
        self.vdim = 3
        self.dpn = 2
    
    def __repr__(self):
        return 'Linear small deformation plane stress elasticity problem'
        
    def NumDofPerNode(self):
        return self.dpn
    
    def SpacialDim(self):
        return self.sdim
        
    def VoigtDim(self):
        return self.vdim
        
    def ComputeInvLumpedMass(self,node,element,dofMap=None):
                
        ndof = dofMap.GID( len(node), 0)
        mv = np.zeros( ndof,basic.FLOAT_TYPE )
        
        etype = None
        mtype = None
        for e in element:
            
            ecoord = node[e.Connectivity(),:]
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                nne = e.NumNodes()
                rho = mtype.matprop['rho'] 
                
            me = np.zeros(nne,basic.FLOAT_TYPE)
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                jac = e.Jacobian(ecoord,qpt)
                N = e.N(qpt)
                me = me + rho*(np.outer(N,N).sum(0))*jac*qwt
           
            for s in xrange(self.sdim):
                sctr = dofMap.Sctr(e.Connectivity(),[s])
                mv[sctr] = mv[sctr] + me
                
        return 1.0/mv
        
        
    def ComputeMass(self,node,element,dofMap=None):
        
        etype = None
        mtype = None
        if ( dofMap == None ):
            dofMap = dmap.FixedDofMap(self.dpn)
            
        #mespace = len(element)*(element[0].NumNodes()*self.dpn)**2
        mdata = basic.DelayedAssm(element,self.dpn)
        
        for e in element:
            ecoord = node[e.Connectivity(),:]
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                nne = e.NumNodes()
                
            if ( not e.Material() == mtype ):
                mtype = e.prop.material
                rho = mtype.matprop['rho'] 
            
            me = np.zeros((nne,nne))
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                jac = e.Jacobian(ecoord,qpt)
                N = e.N(qpt)
                me = me + rho*np.outer(N,N)*jac*qwt
           
            # scatter me into M
            for s in self.SpacialDim():
                sctr = dofMap.Sctr(e.Connectivity(),s)
                mdata.AddLocalMatrix(me,sctr)
                
        return mdata.GetCsrMatrix()
        
    def GetMaterialStiffness(self,e):
        return e.prop.material.TangentStiffness(matl.StressPlaneStress())         
        
    def ComputeStiffness(self,node,element,dofMap=None):
        
        etype = None
        mtype = None
        if ( dofMap == None ):
            dofMap = dmap.FixedDofMap(self.dpn)
        
        #kespace = len(element)*(element[0].NumNodes()*self.dpn)**2
        kdata = basic.DelayedAssm(element,self.dpn)
        
        for e in element:
            ecoord = node[e.Connectivity(),:]
            #pdb.set_trace()
            sctr = dofMap.Sctr(e.Connectivity(),dofMap.LDOFs())
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                kdim = e.NumNodes()*self.dpn
                
            if ( not e.Material() == mtype ):
                mtype = e.prop.material
                C = self.GetMaterialStiffness(e) 

            ke = np.zeros((kdim,kdim))
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                ipm = e.BMat(ecoord,qpt)
                jac = ipm[1]
                B = ipm[0]
                ke = ke + np.dot( np.dot(B.T,C), B )*jac*qwt
           
            # scatter ke into K
            kdata.AddLocalMatrix(ke,sctr)
                
        return kdata.GetCsrMatrix()
         
    
    def ComputeMisesStress(self,sig):
        return np.sqrt( sig[0]**2 + sig[1]**2 + 3*sig[2]**2 - sig[0] * sig[1] )
        
          
    def ComputeStress(self,node,element,disp,dofMap=None):
        
        if ( dofMap == None ):
            dofMap = dmap.FixedDofMap(self.dpn)
            
        ne = len(element)
        strain = np.zeros( (ne,self.vdim), dtype=basic.FLOAT_TYPE )
        stress = np.zeros( (ne,self.vdim), dtype=basic.FLOAT_TYPE )
        svm = np.zeros( ne, dtype=basic.FLOAT_TYPE )
   
        mtype = None
        ei=0
        for e in element:
            ecoord = node[e.Connectivity(),:]
            #sctr = elem.sctr_array( e.Connectivity(), self.dpn )
            sctr = dofMap.Sctr(e.Connectivity(),range(self.dpn))    
            if ( not e.Material() == mtype ):
                mtype = e.prop.material
                C = self.GetMaterialStiffness(e)  
            
            B, jac = e.BMat(ecoord)
            strain[ei,:] = np.dot(B,disp[sctr])
            stress[ei,:] = np.dot(C,strain[ei])
            svm[ei] = self.ComputeMisesStress(stress[ei])
        
            ei = ei + 1
       
        return [ stress, svm, strain ] 

# -------------------------------------------------------------------------        
class ProblemPlaneStrain(ProblemPlaneStress):
    
    def __init__(self):
        self.sdim = 2
        self.vdim = 3
        self.dpn = 2
    
    def __repr__(self):
        return 'Linear small deformation for plane strain problem'
    
    def GetMaterialStiffness(self,e):
        return e.prop.material.TangentStiffness(matl.StressPlaneStrain())         
                                   
# -------------------------------------------------------------------------        
class Problem3D(ProblemPlaneStress):
    
    def __init__(self):
        self.sdim = 3
        self.vdim = 6
        self.dpn = 3
    
    def __repr__(self):
        return 'Linear small deformation for 3D elasticity problem'
    
    def ComputeMisesStress(self,sig):
        return np.sqrt( 0.5*( (sig[0]-sig[1])**2 + (sig[1]-sig[2])**2 + \
                (sig[0]-sig[2])**2 + 6*(sig[3]**2 + sig[4]**2 + sig[5]**2 ) ) )
        
    def GetMaterialStiffness(self,e):
        return e.prop.material.TangentStiffness(matl.Stress3D())         
   
    