import numpy as np
     
   
class PropRod(dict):
    
    def __init__(self,*args,**kwargs):
        dict.__init__(self,*args,**kwargs)
        self.__default_vals__()
         
    def __default_vals__(self):
        param = self.Fields()
        for p in param:
            if ( p not in self ):
                self[p] = param[p]
                
    def Fields(self):
        return { 'PID':None, 'MAT':None, 'A':1.0, 'J':0.0,'C':0.0,'NSM':0.0 } 
                                                                                                                                        
class PropBar(PropRod):
    
    def __init__(self,*args,**kwargs):
        PropRod.__init__(self,*args,**kwargs)
        self.__default_vals__()
        
    def Fields(self):
        return { 'PID':None, 'MAT':None, 
            'A':1.0, 'I1':0.0, 'I2':1.0, 'J':0.0,'NSM':0.0, 
            'SRP':[], 'K1':1.0, 'K2':1.0, 'I12':0.0 }
                                                                                                                                   
class PropShell(PropRod):
    
    def __init__(self,*args,**kwargs):
        PropRod.__init__(self,*args,**kwargs)
        self.__default_vals__()
        
    def Fields(self):
        return { 'PID':None, 'MAT':None, 'T':0.0, 'NSM':0.0,
            'Z1':None, 'Z2':None }
                                                                                                                                   
class PropSolid(PropRod):
    
    def __init__(self,*args,**kwargs):
        PropRod.__init__(self,*args,**kwargs)
        self.__default_vals__()
        
    def Fields(self):
        return { 'PID':None, 'MAT':None }
                   
      
class MaterialPoint(dict):
    
    def __init__(self,*args,**kwargs):
        dict.__init__(self,*args,**kwargs)
 
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
        self.sym = int(sym)
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

                                                     
class LinearElasticMat(dict):
    
    """A class for linear isotropic elastic material"""
    def __init__(self,*args,**kwargs):
        self.mattype = 'linear elastic isotropic'
        dict.__init__(self,*args,**kwargs)
        
    def __repr__(self):
        return self.mattype
        
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
                                                                                                                                         
    def TangentStiffness(self,**kwargs):
        
        try:
            sf = kwargs['sstate']
        except KeyError:
            print 'Error in LinearElasticMat.TangentStiffness: StressState needed'
            return 0
            
        E = self['E']
        G = self['G']
        nu = self['nu'] 
            
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
                    
    def Stress(self,**kwargs):
        strain = kwargs['strain']
        sf = kwargs['sstate']
        C = self.TangentStiffness(sf)
        return np.dot(C,strain)
                                              

               
#------------------------------------------------------------        
#class J2PlasticIsoHard(LinearElasticMat): 
#    """J2 plasticity with isotropic linear hardening"""   
#
#    def __init__(self,*args,**kwargs):
#        self.mattype = 'J2 plastic with linear isotropic hardening'
#        LinearElasticMat.__init__(self,*args,**kwargs)
#                                                                                                               
#    def TangentStiffness(self,**kwargs):
#        
#        try:
#            sf = kwargs['sstate']
#        except KeyError:
#            print 'Error in J2PlasticIsoHard.TangentStiffness: StressState needed'
#            return 0
#            
#    def Stress( self, strain, sf=-1, svar=0 ):
##        import pdb; pdb.set_trace()
#        epsn = strain
#        deps = state.straininc
#        sign = state.stress
#        dsig = self.young*deps
#        if np.abs(sign+dsig) > self.yieldstress :
#            dsig = deps*self.alpha
#            self.yieldstress = self.yieldstress + dsig
#        
#        sig = sign + dsig    
#        state.strain = epsn + deps
#        state.stress = sig  
#        return sig
        

#------------------------------------------------------------ 
#mat = LinearElasticMat(10e6,.28,.028)
#print mat.TangentStiffness()          
#print mat.TangentStiffness(StressPlaneStrain()) 
#print mat.Stress([                                                              
        
     
     
