import numpy as np
                                                                                                                                  
class MaterialState(object):
    """This class abstracts the material state at a point to use in a
    constitutive model to compute stress and other realted material 
    quantities."""
    def __init__(self):
        self.strain=0.0
        self.stress=0.0
        self.straininc=0.0
        
    def __repr__(self):
        return 'strain='+str(self.strain)

#------------------------------------------------------------------------------ 
class Stress3D(object):
    def __init__(self,stress=np.zeros(6,float)):
        self.edim=3   # element dimension
        self.sdim=3   # spacial dim
        self.vdim=6   # voigt dim 
        self.symm=True
        self.SetStress(stress)
        
    def SetStress(self,stress):    
        try:
            self.stress=np.reshape(stress,self.vdim)
        except ValueError:  
            print "Incorrect stress value in Stress __init__"    
        
class StressShell(Stress3D):
    def __init__(self,stress=np.zeros(6,float)):
        self.edim=2
        self.sdim=3
        self.vdim=6     
        self.symm=True
        self.SetStress(stress)
        
class StressRod3D(Stress3D):
    def __init__(self,stress=np.zeros(2,float)):
        self.edim=1
        self.sdim=3
        self.vdim=2   
        self.symm=True
        self.SetStress(stress)
        
class StressPlaneStrain(Stress3D):
    def __init__(self,stress=np.zeros(3,float)):
        self.edim=2
        self.sdim=2
        self.vdim=3   
        self.symm=True
        self.SetStress(stress)
        
class StressPlaneStress(StressPlaneStrain):
    def __init__(self,stress=np.zeros(3,float)):
        self.edim=2
        self.sdim=2
        self.vdim=3    
        self.symm=True
        self.SetStress(stress)
        
class StressAxisymmetric(StressPlaneStrain):
    def __init__(self,stress=np.zeros(4,float)):
        self.edim=2
        self.sdim=2
        self.vdim=4     
        self.symm=True
        self.SetStress(stress)
        
class StressRod2D(Stress3D):
    def __init__(self,stress=np.zeros(1,float)):
        self.edim=1
        self.sdim=2
        self.vdim=1    
        self.symm=True
        self.SetStress(stress)

#------------------------------------------------------------------------------ 
class ConductionMat(object):  
    """A class for linear isotropic conduction material"""
    def __init__(self,kappa,cp,rho):
        self.mattype = 'linear elastic isotropic'
        self.matprop = { 'kappa':kappa, 'cp':cp, 'rho':rho }

#------------------------------------------------------------------------------    
    
class LinearElasticMat(object):
    """A class for linear isotropic elastic material"""
    def __init__(self,E,nu=.3,rho=0):
        self.mattype = 'linear elastic isotropic'
        self.matprop = { 'E':E, 'nu':nu, 'G':E/(1+nu)/2, 'rho':rho }
        
    def __repr__(self):
        return self.mattype
                                                                      
    def TangentStiffness(self,sf=Stress3D(),svar=0):
        
        if ( type(sf) == StressRod2D ):
            C1 = self.matprop['E']
            return np.array([[C1]],float) 
                    
        if ( type(sf) == StressRod3D ):
            C1 = self.matprop['E']
            C2 = self.matprop['G']
            return np.array([[C1,0],[0,C2]],float) 
            
        if ( type(sf) == StressPlaneStress ):
            C1 = self.matprop['E']/(1-(self.matprop['nu'])**2)
            C2 = C1*self.matprop['nu']
            C3 = self.matprop['G']
            return np.array([[C1,C2,0],[C2,C1,0],[0,0,C3]],float)
            
        if ( type(sf) == StressPlaneStrain ):
            C1 = self.matprop['E']/((1+self.matprop['nu'])*(1-2*self.matprop['nu']))
            C2 = C1*self.matprop['nu']
            C1 = (1-self.matprop['nu'])*C1
            C3 = self.matprop['G']
            return np.array([[C1,C2,0],[C2,C1,0],[0,0,C3]],float) 
                    
        if ( type(sf) == StressAxisymmetric ):
            C1 = self.matprop['E']/((1+self.matprop['nu'])*(1-2*self.matprop['nu']))
            C2 = C1*self.matprop['nu']
            C1 = (1-self.matprop['nu'])*C1
            C3 = self.matprop['G']
            return np.array([[C1,C2,C2,0],[C2,C1,C2,0],[C2,C2,C1,0],[0,0,0,C3]],float)
                    
        else:
            C1 = self.matprop['E']/((1+self.matprop['nu'])*(1-2*self.matprop['nu']))
            C2 = C1*self.matprop['nu']
            C1 = (1-self.matprop['nu'])*C1
            C3 = self.matprop['G']
            return np.array([[C1,C2,C2,0,0,0],[C2,C1,C2,0,0,0],\
                    [C2,C2,C1,0,0,0],[0,0,0,C3,0,0],[0,0,0,0,C3,0],\
                    [0,0,0,0,0,C3]],float) 
                    
    def Stress( self, strain ):
        C = self.TangentStiffness(strain)
        #if ( type(strain) != numpy.ndarray ):
        #    return np.dot(C,array(strain).reshape((6,1))) 
        #return np.dot(C,strain.reshape((6,1)))    
               
#------------------------------------------------------------        
class J2PlasticIsoHard(LinearElasticMat): 
    """J2 plasticity with isotropic linear hardening"""   
    def __init__(self,E,nu=.3,rho=0,sy=0,alpha=0):
        self.young = E
        self.poisson = nu
        self.density = rho
        self.yieldstress = sy
        self.alpha = alpha 
        
    def Stress( self, strain, sf=-1, svar=0 ):
#        import pdb; pdb.set_trace()
        epsn = strain
        deps = state.straininc
        sign = state.stress
        dsig = self.young*deps
        if np.abs(sign+dsig) > self.yieldstress :
            dsig = deps*self.alpha
            self.yieldstress = self.yieldstress + dsig
        
        sig = sign + dsig    
        state.strain = epsn + deps
        state.stress = sig  
        return sig
        

#------------------------------------------------------------ 
#mat = LinearElasticMat(10e6,.28,.028)
#print mat.TangentStiffness()          
#print mat.TangentStiffness(StressPlaneStrain()) 
#print mat.Stress([                                                              