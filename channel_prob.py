import numpy as np
import basic
import gmsh
import dofmap as dmap
import material as matl
import prop
 
                                
# ------------------------------------------------------------------------- 
class CoolingChannelProblem(object):
    
    def __init__(self,dimension=3,dofpn=1):
        self.sdim = dimension
        self.dpn = dofpn     
        
    def ComputeStiffness(self,elements,node,dofMap):
        
        kdata = basic.DelayedAssm(block,self.dpn)
        
        etype = None
        mtype = None
        
        for block in elements:
            
            for e in block:
                
                if ( not etype==e.Type() ):
                    etype = e.Type()
                    qr = e.QuadratureRule(2*(e.Order()))
                    nne = e.NumNodes()
                  
                if ( not e.Material() == mtype ):
                    mtype = e.prop.material
                    C = self.GetMaterialStiffness(e) 
                
                ecoord = node[e.Connectivity(),:]
                sctr = dofMap.Sctr(e.Connectivity(),dofMap.LDOFs())
            
                ke = np.zeros((nne,nne))
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
                                    
                                                          
# -----------------------------------
gmshfile = gmsh.GmshInput('channel.msh')
channelPID = 163
fluidPID = 162

inletPID = 164
heatFluxPID = 154

gmshfile.AddPhysicalIDs([channelPID,fluidPID])

gmshfile.AddSideSetIDs( inletPID )
gmshfile.AddNodeSetIDs( heatFluxPID )

elements = gmshfile.ReadElements()  
[nodes,nids] = gmshfile.ReadNodes()    
    
nodesets = gmshfile.ReadNodeSets()
sidesets = gmshfile.ReadSideSets()

dofmap = dmap.VariDofMap(nids)
dofmap.ActivateDofs(elements,1)
dofmap.FixDofs()
    
wallMat = matl.LinearElasticMat(E=10e6,nu=.3,rho=.098,kappa=236,cp=10,alpha=2e-6)
wallProp = prop.PropSolid(Matl=wallMat)    
    
