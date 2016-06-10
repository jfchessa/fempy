import numpy as np
import basic
import element as elem  
import dofmap as dmap
import linear_problems as prob
import bcs
import sets
import pdb
              
# -------------------------------------------------------------------------   
class NewtonianFluid(object):
    
    def __init__(self,rho=1.003,mu=.0024):
        self.matprop = {'rho':rho, 'mu':mu}
        
    def ShearStress(self,D):
        Dh=D.trace()/len(D)
        Ddev = D
        for s in xrange(len(D)):
            Ddev[s,s] = Ddev[s,s] - Dh
        return 2*self.matprop['rho']*Ddev
    
class CBSMethod(prob.Problem3D):
    
    def __init__( self, sdim=3 ):
        self.sdim = sdim
        self.dpn = self.sdim+1
        self.grav= np.zeros( self.sdim, basic.FLOAT_TYPE )
    
    def GetMaterialStiffness(self,e):
        mu = e.prop.material.matprop['mu'] 
        if ( self.sdim==3 ):
            C = np.zeros( (6,6), basic.FLOAT_TYPE )
            C[:3,:3] = -2*mu/3
            C[0,0] = 4*mu/3
            C[1,1] = C[0,0]
            C[2,2] = C[0,0]
            C[3,3] = mu
            C[4,4] = mu
            C[5,5] = mu
            return C
            
        elif ( self.sdim==2 ):
            C = np.zeros( (3,3), basic.FLOAT_TYPE )
            C[:2,:2] = -2*mu/3
            C[0,0] = 4*mu/3
            C[1,1] = C[0,0]
            C[2,2] = mu
            return C
            
        else: # sdim ==1
            return np.array( [2*mu], basic.FLOAT_TYPE )
                                 
    def ComputeInvMu(self,node,element,dofmap):       
        ndof = dofmap.GID( len(node), 0)
        mv = np.zeros( ndof,basic.FLOAT_TYPE )
        
        etype = None
        for e in element:
            
            ecoord = node[e.Connectivity(),:]
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                nne = e.NumNodes()
                
            me = np.zeros(nne,basic.FLOAT_TYPE)
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                jac = e.Jacobian(ecoord,qpt)
                N = e.N(qpt)
                me = me + (np.outer(N,N).sum(0))*jac*qwt
           
            for s in xrange(self.sdim):
                sctr = dofmap.Sctr(e.Connectivity(),[s])
                mv[sctr] = mv[sctr] + me
                
        return 1.0/mv

            
    def Step1RHS(self,node,element,dofmap,U,dt):
        ndof = dofmap.GID( len(node), 0)
        rhs = np.zeros( ndof,basic.FLOAT_TYPE )
        
        etype = None
        mtype = None
        for e in element:
            
            ecoord = node[e.Connectivity(),:]
            sctr = dofmap.Sctr(e.Connectivity())
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                nne = e.NumNodes()
            if ( not mtype==e.Material() ):
                mtype = e.Material()
                rhoe = mtype.matprop['rho']
                
            re = np.zeros(self.sdim*nne,basic.FLOAT_TYPE)
            Ui = U[sctr].reshape(nne,self.sdim)
            
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                
                N = e.N(qpt)
                dNdx, jac = e.dNdx(ecoord,qpt)
                
                Upt = np.dot(N,Ui)
                gradU = np.dot(dNdx.T,Ui)
                divU = gradU.trace()
                upt = Upt/rhoe
                gradu = gradU/rhoe
                divu = gradu.trace()
                divuU = divu*Upt + upt*divU
                tau = e.Material().ShearStress(gradu)
                Btau = np.dot(dNdx,tau.T)
                    
                for i in xrange(self.sdim):
                    re[i::self.sdim] = re[i::self.sdim] + ( -N*divuU[i] - Btau[:,i] + N*rhoe*self.grav[i] + \
                        0.5*dt*( (divu*N+np.dot(dNdx,upt))*(-divuU[i] + rhoe*self.grav[i]) ))*jac*qwt

                rhs[sctr] = rhs[sctr] + re
                
        return rhs
    
    def Step2Operators(self,node,elemp,dofmap,elemu=None):
         
        if elemu==None:
            elemu = elemp

        mpdata = basic.DelayedAssm(elemp,1)
        kpdata = basic.DelayedAssm(elemp,1)
        gpdata = basic.DelayedAssm(elemp,self.sdim)

        etype = None
        ee=0
        for e in elemp:
            
            ecoord = node[e.Connectivity(),:]
            if ( not etype==e.Type() ):
                etype = e.Type()
                qr = e.QuadratureRule(2*(e.Order()))
                nne = e.NumNodes()
                
            me = np.zeros((nne,nne),basic.FLOAT_TYPE)
            ke = np.zeros((nne,nne),basic.FLOAT_TYPE)
            ge = np.zeros((nne,self.sdim*nne),basic.FLOAT_TYPE)
            
            for q in xrange(qr[1].size):
                qpt = qr[0][q]
                qwt = qr[1][q]
                
                N = e.N(qpt)
                Nu = elemu[ee].N(qpt)
                dNdx, jac = e.dNdx(ecoord,qpt)
                c=1.0  # ????^-1QUE????
               
                me = me + (1/c**2)*(np.outer(N,N))*jac*qwt
                ke = ke + np.dot(dNdx,dNdx.T)*jac*qwt
                for i in xrange(self.sdim):
                    ge[:,i::self.sdim] = ge[:,i::self.sdim] + (np.outer(dNdx[:,i],Nu))*jac*qwt
           
            sctru = dofmap.Sctr(e.Connectivity())
            sctrp = e.Connectivity()
            
            mpdata.AddLocalMatrix(me,sctrp)
            kpdata.AddLocalMatrix(ke,sctrp)
            gpdata.AddLocalMatrix(ge,sctrp,sctru)
            ee=ee+1
            
        return [ mpdata.GetCsrMatrix(), kpdata.GetCsrMatrix(), gpdata.GetCsrMatrix() ]
          
                           
        def Step2Influx(self,node,inflowbc):
            return []
            
class VBC(object):
    def __init__(self,felem,n,idof=None,ival=None):
        self.felem = felem
        self.normal = n
        if ( idof == None ):
            self.idof = np.array([0,1,2],basic.INDX_TYPE)
        else:
            self.idof = idof
        if ( ival == None ):
            self.ival = 0*idof
        else:
            self.ival = ival
                
class VelocityBCs(object):
    
    def __init__(self):
       self.vbcs =[]
       self.spcs=[]
       
       
    def AddBC(self,element,n,idof=None,ival=None):
        ee=0
        nid = sets.Set()
        for e in element:
            self.vbcs.append(VBC(e,n[ee],idof,ival))
            ee = ee + 1
            for I in e.conn:
                for i in xrange(len(idof)):
                    nid.add( (I,idof[i],ival[i]) )
                
        for s in nid:
            self.spcs.append(bcs.SPC(s[0],s[1],s[2]))
        
        
    def SetValue(self,U,dofmap):
                
        for s in self.spcs:
            Ii = dofmap.GID(s.nid,s.ldof)
            U[Ii] = s.dval
        return U
                  
    def GetPPFlux(self,node,element,dofmap):       
        nn = len(node)  #  need to rework this whole thing
        ndof = dofmap.GID(nn,0)
        fext = np.zeros( ndof, basic.FLOAT_TYPE )
                
        for v in self.vbcs:
            
            dpn = len(f.trac)
            fe = np.zeros( f.face.NumNodes()*dpn, basic.FLOAT_TYPE )
            ldof = np.array( range(dpn), basic.INDX_TYPE )
            sctr = dofmap.Sctr(f.face.Connectivity(),ldof)
            
            # compute fe
            etype=None
            e = f.face
            nvect = v.normal  
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
                fe = fe + N*v.ival*nvect[v.idof]*jac*qwt
                    
            fext[sctr] = fext[sctr] + fe
            
        return fext
          
                                                                                                                  
#---------------------------------------------
import prop
import meshing as mesh

# define the finite element mesh
pwidth=20
pheight=5
mat0 = NewtonianFluid( )
prop0 = prop.PlaneStress( mat0, 1.0 )
fegrid = mesh.MeshQuad4( np.array([[0,0],[pwidth,0],[pwidth,pheight],[0,pheight]]), 5, 2, prop0 )
dofmap = dmap.FixedDofMap(2)

# define the velocity boundary conditions
inflow = bcs.NaturalBCs()
inflow.AddTraction( fegrid.EdgeElem(3), [10.0, 0.0] )

inflowspcs = bcs.EssentialBCs()
inflowspcs.AddSpcs( fegrid.EdgeNIDs(3), [0], [10.0], dofmap=dofmap )  # fix left edge in x
loads = bcs.NaturalBCs()  # if there are external wall tractions on the flow
loads.AddTraction( fegrid.EdgeElem(0), [10.0, 0.0] )

# set the initial conditions
formulation = CBSMethod(sdim=2) 
dt=.1
U=np.zeros(20,float)
U[::2]=fegrid.node[:,0]

# solve a single time increment
# step 1
mlinv = formulation.ComputeInvMu(fegrid.node,fegrid.element,dofmap) 
re = formulation.Step1RHS(fegrid.node,fegrid.element,dofmap,U,dt)
fext = loads.Fext(fegrid.node,dofmap)  
dU = dt*mlinv*(re+fext)
U = U + dU
#spcs.SetValues(U)     

# step 2
mp, kp, gp = formulation.Step2Operators(fegrid.node,fegrid.element,dofmap)
#ppsolver = basic.FeSolver()
norms = fegrid.EdgeNormals(0)

# step 3

#solve 

# step 3     