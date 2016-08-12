import numpy as np
              
 
class SS3DStress(object):
    def __init__(self):
         self.vdim=6
         self.dim=3
         self.ndi=3
         self.nshr=3
          
class SSPlaneStress(object):
    def __init__(self):
         self.vdim=3
         self.dim=2
         self.ndi=2
         self.nshr=1
                              
class SSPlaneStrain(object):
    def __init__(self):
         self.vdim=3
         self.dim=2
         self.ndi=2
         self.nshr=1
                                                
class SSAxisymmetric(object):
    def __init__(self):
         self.vdim=3
         self.dim=2
         self.ndi=2
         self.nshr=1
                          
class SSRod2D(object):
    def __init__(self):
         self.vdim=1
         self.dim=2
         self.ndi=1
         self.nshr=0
                                    
class SSRod3D(object):
    def __init__(self):
         self.vdim=2
         self.dim=2
         self.ndi=1
         self.nshr=1

class MaterialPoint(dict):
    
    def __init__(self,*args,**kwargs):
        dict.__init__(self,*args,**kwargs)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
class Tensor(np.ndarray):
 
    def __new__( subtype, rank, ndim, sym=False, buff=None ): 
        
        if ( rank==0 ):
            shp = 1
        elif ( sym ):
            if ( rank == 2 ):
                shp = ( ndim**rank - ndim )/2 + ndim
            elif ( rank == 4 ):
                r=rank/2
                shp = ( ( ndim**r - ndim )/2 + ndim, )*r
            else:
                print 'Error, rank of symetric tensor must be either 2 or 4' 
                sym =False    
                  
        if ( not sym ):
            shp = (ndim,)*rank
 
        try:
            obj = np.ndarray.__new__(subtype, shp, float, np.array(buff,float) )
        except:
            buff = None
            obj = np.ndarray.__new__(subtype, shp, float )
            
        return obj
        
    def __init__( self, rank, ndim, sym=False, buff=None  ):
        self.sym = sym
        self.rank = int(rank)
        self.dim = int(ndim)
        if ( buff==None ):
            self.fill(0)

    def Rank(self):
        return self.rank;
        
    def Dim(self):
        return self.dim   
         
    def Set( self, buff ):
        self = Tensor( self.Dim(), buff )
                              
    #def Full( self ):
    #    
    #    return self.reshape((1,9))[ 0, [0,4,8,5,2,1] ] 
   
class SR2Tensor1D(Tensor):
              
    def __new__( subtype, buff=None ): 
        obj = Tensor.__new__(subtype,2,1,True,buff)   
        return obj
                                                    
    def __init__( self, buff=None  ):
        Tensor.__init__(self,2,1,True,buff) 
        
    def Hydrostatic(self):
        return self
    
    def Deviatoric(self):
        return SR2Tensor1D([0.0])
        
    def Mises(self):
        return abs(self)        
        
class SR2Tensor2D(Tensor):
              
    def __new__( subtype, buff=None ): 
        obj = Tensor.__new__(subtype,2,2,True,buff)   
        return obj
                                                    
    def __init__( self, buff=None  ):
        Tensor.__init__(self,2,2,True,buff) 
        
    def Hydrostatic(self):
        shyd = 0.5*(self[0] + self[1])
        return SR2Tensor2D( [shyd, shyd, 0.0] )
    
    def Deviatoric(self):
        return self - self.Hydrostatic()
        
    def Mises(self):
        vm = (self[0]-self[1])**2 + self[0]**2 + self[1]**2 + 6*self[2]**2
        return np.sqrt(0.5*vm)
                
class SR2Tensor3D(Tensor):
              
    def __new__( subtype, buff=None ): 
        obj = Tensor.__new__(subtype,2,3,True,buff)   
        return obj
                                                    
    def __init__( self, buff=None  ):
        Tensor.__init__(self,2,3,True,buff)   
    
    def Hydrostatic(self):
        shyd = 0.333333333333333333333*(self[0]+self[1]+self[2])
        return SR2Tensor3D( [shyd, shyd, shyd, 0.0, 0.0, 0.0] )
    
    def Deviatoric(self):
        return self - self.Hydrostatic()
        
    def Mises(self):
        vm = (self[0]-self[1])**2 + (self[1]-self[2])**2 + (self[2]-self[0])**2
        vm += 6*(self[3]**2 + self[4]**2 + self[5]**2)
        return np.sqrt(0.5*vm)
                              

def SR2Tensor(dim,buff=None):
    if ( dim == 1 ):
        return SR2Tensor1D(buff)
    elif ( dim == 2 ):
        return SR2Tensor2D(buff)
    else:
        return SR2Tensor3D(buff)
        
                                                                                                                                                                                                                                                                                                                                    
class LinearElasticMat(dict):
    
    """A class for linear isotropic elastic material"""
    def __init__(self,*args,**kwargs):
        self.mattype = 'linear elastic isotropic'
        dict.__init__(self,*args,**kwargs)
        self.CheckProperties()
        
    def __repr__(self):
        return self.mattype
     
    def InitializeMaterialPoint(self,mp):
        sf = mp['StressState']
        sd = sf.ndi + sf.nshr
        mp['Stress'] = SR2Tensor(sd)
           
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
            
        if ( type(sf) == SSRod2D ):
            return np.array([[E]],float)
                    
        if ( type(sf) == SSRod3D ):
            return np.array([[E,0],[0,G]],float)  
            
        if ( type(sf) == SSPlaneStress ):
            C1 = E/(1-nu**2)
            C2 = C1*nu
            C3 = G
            return np.array([[C1,C2,0],[C2,C1,0],[0,0,C3]],float)
            
        if ( type(sf) == SSPlaneStrain ):
            C1 = E/((1+nu)*(1-2*nu))
            C2 = C1*nu
            C1 = (1-nu)*C1
            C3 = G
            return np.array([[C1,C2,0],[C2,C1,0],[0,0,C3]],float) 
                    
        if ( type(sf) == SSAxisymmetric ):
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
        mp['Strain']    = SR2Tensor(sd)
        mp['Stress']    = SR2Tensor(sd)
        mp['StrainInc'] = SR2Tensor(sd)
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
mat0 = LinearElasticMat(E=10e6,nu=.28)
mp = MaterialPoint( StressState=SSRod2D() )
C = mat0.TangentStiffness( mp )    

mat1 = J2PlasticIsoHard(E=10e6,nu=.28,Sy=10000,Etan=5e5) 
mat1.InitializeMaterialPoint(mp)  
                                         
mp['StrainInc'] =  SR2Tensor1D([.0001])    
sig=[]
eps=[]
for i in xrange(25):
    sig.append( mat1.Stress(mp).Mises() )
    eps.append( mp['Strain'].Mises() )  
     