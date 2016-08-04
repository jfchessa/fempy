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
        
        kdata = basic.DelayedAssm(elements,self.dpn)
        
        etype = None
        mtype = None
        
        for e in elements:
            
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                nne = e.NumNodes()
                
            if ( not e.Material() == mtype ):
                mtype = e.prop['Matl']
                k = mtype['kappa'] 
            
            ecoord = node[e.Connectivity(),:]
            sctr = dofMap.Sctr(e.Connectivity(),[0])
        
            ke = np.zeros((nne,nne))
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                [dNdx,jac] = e.dNdx(ecoord,qpt)
                ke = ke + np.dot( dNdx, dNdx.T )*k*jac*qwt
        
            # scatter ke into K
            kdata.AddLocalMatrix(ke,sctr)
                
        return kdata.GetCsrMatrix()     
                                    
                                                          
# ---------------------------------------------------------------------------
gmshfile = gmsh.GmshInput('channel.msh')

channelPID = 163
fluidPID = 162

inletPID = 164
heatFluxPID = 154

wallMat = matl.LinearElasticMat(E=10e6,nu=.3,rho=.098,kappa=236,cp=10,alpha=2e-6)
fluidMat = matl.LinearElasticMat(E=10e6,nu=.3,rho=.098,kappa=236,cp=10,alpha=2e-6)
wallProp = prop.PropSolid(Matl=wallMat,Region='Solid')
fluidProp = prop.PropSolid(Matl=fluidMat,Region='Fluid',Avel=1.0)
props = { channelPID:wallProp, fluidPID:fluidProp }  

#gmshfile.AddPhysicalIDs([channelPID,fluidPID])
gmshfile.AddPhysicalIDs([channelPID])

gmshfile.AddSideSetIDs( inletPID )
gmshfile.AddNodeSetIDs( heatFluxPID )

elements = gmshfile.ReadElements(props)  
[nodes,nids] = gmshfile.ReadNodes()   
elements = gmshfile.RenumberConnectivity(nids,elements) 
    
nodesets = gmshfile.ReadNodeSets()
sidesets = gmshfile.ReadSideSets()

dofmap = dmap.VariDofMap(nids)
dofmap.ActivateDofs(elements,1)
dofmap.FixDofs()
    
problem = CoolingChannelProblem()
K = problem.ComputeStiffness(elements,nodes,dofmap)
