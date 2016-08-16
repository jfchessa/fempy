import numpy as np
import mat_basic as mb
                                                                                                                                                                                                                                                                                                       
class LinearElasticMat(dict):
    
    """A class for linear isotropic elastic material 
     LinearElasticMat(E=,nu=,rho=...)
    """
    def __init__(self,*args,**kwargs):
        self.mattype = 'linear elastic isotropic'
        dict.__init__(self,*args,**kwargs)
        self.CheckProperties()
        
    def __repr__(self):
        return self.mattype
     
    def InitializeMaterialPoint(self,mp):
        sf = mp['StressState']
        sd = sf.ndi + sf.nshr
        mp['Stress'] = mb.SR2Tensor(sd)
           
    def CheckProperties(self):
        if ( ('E' in self) and ('nu' in self) ):
            self['G'] = 2.0*self['E']/(1+self['nu'])
            
        elif ( ('E' in self) and ('G' in self) ):
            self['nu'] = 2*self['E']/self['G'] - 1.0
            
        elif ( ('nu' in self) and ('G' in self) ):
            self['E'] = 0.5*self['G']*(1.0+self['nu'])
            
        else:
            return True
        
        return False
                                                                                                                                         
    def TangentStiffness(self,mp):
            
        E = self['E']
        G = self['G']
        nu = self['nu'] 
        
        sf = mp['StressState']    
            
        if ( type(sf) == mb.SSRod2D ):
            return np.array([[E]],float)
                    
        if ( type(sf) == mb.SSRod3D ):
            return np.array([[E,0],[0,G]],float)  
            
        if ( type(sf) == mb.SSPlaneStress ):
            C1 = E/(1-nu**2)
            C2 = C1*nu
            C3 = G
            return np.array([[C1,C2,0],[C2,C1,0],[0,0,C3]],float)
            
        if ( type(sf) == mb.SSPlaneStrain ):
            C1 = E/((1+nu)*(1-2*nu))
            C2 = C1*nu
            C1 = (1-nu)*C1
            C3 = G
            return np.array([[C1,C2,0],[C2,C1,0],[0,0,C3]],float) 
                    
        if ( type(sf) == mb.SSAxisymmetric ):
            C1 = E/((1+nu)*(1-2*nu))
            C2 = C1*nu
            C1 = (1-nu)*C1
            C3 = G
            return np.array([[C1,C2,C2,0],[C2,C1,C2,0],[C2,C2,C1,0],[0,0,0,C3]],float)
                    
        else:
            C1 = E/((1+nu)*(1-2*nu))
            C2 = C1*nu
            C1 = (1-nu)*C1
            C3 = G
            return np.array([[C1,C2,C2,0,0,0],[C2,C1,C2,0,0,0],\
                    [C2,C2,C1,0,0,0],[0,0,0,C3,0,0],[0,0,0,0,C3,0],\
                    [0,0,0,0,0,C3]],float)   
                    
    def Stress(self,mp):
        sf = mp['StressState']    
        strain = mp['Strain']
        C = self.TangentStiffness(sf)
        return np.dot(C,strain)

          
                                                
class J2PlasticIsoHard(LinearElasticMat): 
    """J2 plasticity with isotropic linear hardening"""   

    def __init__(self,*args,**kwargs):
        """ needs E, Sy and Etan defined """
        LinearElasticMat.__init__(self,*args,**kwargs)
        self.mattype = 'J2 plastic with linear isotropic hardening'
        
    def InitializeMaterialPoint(self,mp):
        sf = mp['StressState']
        sd = sf.ndi + sf.nshr
        mp['Strain']    = mb.SR2Tensor(sd)
        mp['Stress']    = mb.SR2Tensor(sd)
        mp['StrainInc'] = mb.SR2Tensor(sd)
        mp['Sy'] = self['Sy']
                                                                                                                       
    def ElasticTanStiffness(self,mp):
        return LinearElasticMat.TangentStiffness(self,mp) 
                                                                                           
    def InelasticTanStiffness(self,mp):
        C = LinearElasticMat.TangentStiffness(self,mp)
        return (self['Etan']/self['E'])*C
                                                                                            
    def TangentStiffness(self,mp):
        
        svm = mp['Stress'].Mises()
        C = LinearElasticMat.TangentStiffness(self,mp)
        if ( svm < mp['Sy'] ):
            return C
        else:
            return (self['Etan']/self['E'])*C
            
    def Stress( self, mp ):
        """Requires Stress, Strain, StrainInc, and Sy in the MaterialPoint"""
#        import pdb; pdb.set_trace()
        epsn = mp['Strain']
        deps = mp['StrainInc']
        sign = mp['Stress']
     
        yieldstress = mp['Sy']
        
        dsig = self.ElasticTanStiffness(mp)*deps
        strial = sign + dsig
        if ( strial.Mises() >= yieldstress ):
            dsig = self.InelasticTanStiffness(mp)*deps
            mp['Sy'] = mp['Sy'] + dsig.Mises()
        
        sig = sign + dsig    
        mp['Strain'] = epsn + deps
        mp['Stress'] = sig  
        return sig
        

#------------------------------------------------------------ 
#mat0 = LinearElasticMat(E=10e6,nu=.28)
#mp = mb.MaterialPoint( StressState=mb.SSRod2D() )
#C = mat0.TangentStiffness( mp )    
#
#mat1 = J2PlasticIsoHard(E=10e6,nu=.28,Sy=10000,Etan=5e5) 
#mat1.InitializeMaterialPoint(mp)  
#                                         
#mp['StrainInc'] =  mb.SR2Tensor1D([.0001])    
#sig=[]
#eps=[]
#for i in xrange(25):
#    sig.append( mat1.Stress(mp).Mises() )
#    eps.append( mp['Strain'].Mises() )  
#     