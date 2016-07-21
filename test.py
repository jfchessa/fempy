import numpy as np
     
                                                                                                          
class Tensor(np.ndarray):
    
    def __new__( subtype, ndim, buff=0 ):  
        try:
            obj = np.ndarray.__new__(subtype, (ndim,ndim), float, np.array(buff,float) )
        except:
            obj = np.ndarray.__new__(subtype, (ndim,ndim), float, np.zeros((ndim,ndim)) )
            
        return obj
        
    def __init__( self, ndim, buff=0  ):
        self.sym =False
        
    def Rank(self):
        return 2;
        
    def Dim(self):
        return int(self.shape[0])   
         
    def Set( self, buff ):
        self = Tensor( self.Dim(), buff )
        
class Tensor3D(Tensor):
    
    def __new__( subtype, buff=0 ):   
        obj = Tensor.__new__( subtype, 3, buff )
        return obj
        
    def __init__( self,  buff=0  ):
        self.sym =False
                   
    def voigt(self):
        return self.reshape((1,9))[ 0, [0,4,8,5,2,1] ]
                                                                                                
class Tensor2D(Tensor):
    
    def __new__( subtype, buff=0 ):  
        obj = Tensor.__new__( subtype, 2, buff )
        return obj
          
    def __init__( self, buff=0  ):
        self.sym =False  
                   
    def voigt(self):
        return self.reshape((1,6))[ 0, [0,3,1] ]
                                                                                                 
class Tensor1D(Tensor):
    
    def __new__( subtype, buff=0 ):   
        obj = Tensor.__new__( subtype, 1, buff )
        return obj
        
    def __init__( self, buff=0 ):
        self.sym =False
        
    def voigt(self):
        return self.reshape((1,1))[ 0, 0 ]

#-------------------------------------------------------------------------------    
class Stress3D(Tensor3D):
    
    def __init__( self, buff=0 ):
        self.vdim=6   
        
class StressPlaneStress(Tensor2D):
    
    def __init__( self, buff=0 ):
        self.vdim=3   

        
     
     
