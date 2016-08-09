import numpy as np
import basic
import gmsh
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
        
        for eid, e in elements.iteritems():
            
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                nne = e.NumNodes()
                
            if ( not e.Material() == mtype ):
                mtype = e.prop['Matl']
                k = mtype['kappa'] 
            
            ecoord = node[ e.Connectivity() ]
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
heatFluxPID = 165

wallMat = matl.LinearElasticMat(E=10e6,nu=.3,rho=.098,kappa=236,cp=10,alpha=2e-6)
fluidMat = matl.LinearElasticMat(E=10e6,nu=.3,rho=.098,kappa=236,cp=10,alpha=2e-6)
wallProp = prop.PropSolid(Matl=wallMat,Region='Solid')
fluidProp = prop.PropSolid(Matl=fluidMat,Region='Fluid',Avel=1.0)
props = { channelPID:wallProp, fluidPID:fluidProp }  

#gmshfile.AddPhysicalIDs([channelPID,fluidPID])
gmshfile.AddPhysicalIDs([channelPID])

gmshfile.AddSideSetIDs( heatFluxPID )
gmshfile.AddNodeSetIDs( inletPID )

nodes = gmshfile.ReadNodes()   
elements = gmshfile.ReadElements(props)   
    
nodesets = gmshfile.ReadNodeSets()
sidesets = gmshfile.ReadSideSets()

dofmap = basic.DofMap()
dofmap.ActivateDofs(elements,1)
dofmap.Renumber()

problem = CoolingChannelProblem()
K = problem.ComputeStiffness(elements,nodes,dofmap)

spcs = basic.EssentialBCs()
spcs.AddPointValues(nodesets[inletPID],[0],[200.0])

force = basic.NaturalBCs()
force.AddFaceTraction(nodes,sidesets[heatFluxPID],[0],[10.0])

fext = np.zeros( dofmap.NumDof(), basic.FLOAT_TYPE )
force.AddRHS(dofmap,fext)

[ifix,ival] = spcs.GetIFIX(dofmap)
[d,freac] = basic.fesolve(K,fext,ifix,ival)
