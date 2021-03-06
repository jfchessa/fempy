import numpy as np
import fempy as fp
import fempy.meshing as msh
import fempy.material as mat
import fempy.io as io

#-------------------------------------------------------------------------------
class sol101(object):
    
    def __init__(self):
        self.dofpn=3
        self.sdim=3
        
    
    def GetCMat(self,e):
        
        mp = mat.MaterialPoint()
        mp['StressState'] = mat.SS3DStress()
        return e.prop['Matl'].TangentStiffness(mp)
        
    def ComputeStiffnessMat(self,node,element,dofmap):
        
        etype = None
        mtype = None
        
        kdata = fp.DelayedAssm(element,self.dofpn)
        
        for e in element.values():
            ecoord = node[e.Connectivity()]

            sctr = dofmap.Sctr(e.Connectivity(),range(self.dofpn))
            if ( not etype==type(e) ):
                etype = type(e)
                qr = e.QuadratureRule(2*(e.Order()))
                kdim = e.NumNodes()*self.dofpn
               
            if ( not type(e.prop['Matl']) == mtype ):
                mtype = type(e.prop['Matl'])
                C = self.GetCMat(e)

            ke = np.zeros((kdim,kdim))
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                ipm = e.BMat(ecoord,qpt)
                jac = ipm[1]
                B = ipm[0]
                ke = ke + np.dot( np.dot(B.T,C), B )*(jac*qwt)
           
            # scatter ke into K
            kdata.AddLocalMatrix(ke,sctr)
                
        return kdata.GetCsrMatrix()

    def PostProcess(self,filename,node,element,dofmap,dd):
        
        mtype = None
        ne = len(element)
        stress = np.zeros( (ne,9), float )
        svm = np.zeros(ne,float)
        ee=0
        for eitr in element.iteritems():
            
            e = eitr[1]
            
            ecoord = node[e.Connectivity()]

            sctr = dofmap.Sctr(e.Connectivity(),range(self.dofpn))
               
            if ( not type(e.prop['Matl']) == mtype ):
                mtype = type(e.prop['Matl'])
                C = self.GetCMat(e)

            xi = e.QuadratureRule(1)[0][0] 
            ipm = e.BMat(ecoord,xi)
            B = ipm[0]
            estrain = np.dot(B,dd[sctr])
            estress = np.dot(C,estrain)
            if len(estress)==3:
                stress[ee,[0,1,3]]=estress
            else:
                stress[ee,:len(estress)]=estress
            svm[ee] = np.sqrt( 0.5*(  (stress[ee,0]-stress[ee,1])**2 + \
                (stress[ee,1]-stress[ee,2])**2 + (stress[ee,2]-stress[ee,0])**2 + \
                6*(stress[ee,3]**2 +stress[ee,4]**2+stress[ee,5]**2)  ) )
            ee += 1
                
        dataout = io.FeaData(mesh.node,mesh.element)
        dataout.SetDisplacement(dd,dofmap)
        dataout.SetStress(stress)
        dataout.SetEffectiveStress(svm)
        dataout.WriteVtkFile(filename)

#-------------------------------------------------------------------------------
L=10
W=5
H=2

matl = mat.LinearElasticMat(E=10.0e6,nu=.3)
prop = fp.PropSolid(Matl=matl)

corners = np.array([[0,0,0],[L,0,0],[L,W,0],[0,W,0],[0,0,H],[L,0,H],[L,W,H],[0,W,H]])
mesh = msh.MeshHexa8( corners, 51, 11, 5, prop )

dofmap = fp.DofMap()
dofmap.ActivateDofs(mesh.element,[0,1,2])

spcs = fp.EssentialBCs()
spcs.AddPointValues(mesh.FaceNIDs(1),0)
spcs.AddPointValues(mesh.FaceNIDs(3),1)
spcs.AddPointValues(mesh.FaceNIDs(5),2)

load = fp.NaturalBCs()
tface = mesh.FaceElem(0)
load.AddFaceTraction(mesh.node,mesh.FaceElem(0),[0],[20.0])

f = np.zeros(dofmap.NumDof(),float)
load.AddRHS(dofmap,f)

problem = sol101()
K = problem.ComputeStiffnessMat(mesh.node,mesh.element,dofmap)

ifix = spcs.GetIFIX(dofmap)
[d,freac]=fp.fesolve(K,f,ifix[0])

problem.PostProcess("ex1.vtu",mesh.node,mesh.element,dofmap,d)