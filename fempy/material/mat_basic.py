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
                          
class SSRod1D(object):
    def __init__(self):
         self.vdim=1
         self.dim=1
         self.ndi=1
         self.nshr=0
  
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
        
           