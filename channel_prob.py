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
        
    def ComputeStiffness(self,elements,dofmap):
        
        etype = None
        mtype = None
        for block in elements:
            
            kdata = basic.DelayedAssm(block,self.dpn)
            for e in block:
                
                if ( not etype==e.Type() ):
                    etype = e.Type()
                    qr = e.QuadratureRule(2*(e.Order()))
                    nne = e.NumNodes()
                    
                                    
                                                          
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
dofmap.ActivateDofs(elements[channelPID],1)
dofmap.ActivateDofs(elements[fluidPID],1)
dofmap.FixDofs()
    
wallMat = matl.LinearElasticMat(E=10e6,nu=.3,rho=.098,kappa=236,cp=10,alpha=2e-6)
wallProp = prop.PropSolid(Matl=wallMat)    
    
