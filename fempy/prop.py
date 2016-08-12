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
        return { 'PID':None, 'Matl':None, 'A':1.0, 'J':0.0,'C':0.0,'NSM':0.0 } 
                                                                                                                                        
class PropBar(PropRod):
    
    def __init__(self,*args,**kwargs):
        PropRod.__init__(self,*args,**kwargs)
        self.__default_vals__()
        
    def Fields(self):
        return { 'PID':None, 'Matl':None, 
            'A':1.0, 'I1':0.0, 'I2':1.0, 'J':0.0,'NSM':0.0, 
            'SRP':[], 'K1':1.0, 'K2':1.0, 'I12':0.0 }
                                                                                                                                      
class PropShell(PropRod):
    
    def __init__(self,*args,**kwargs):
        PropRod.__init__(self,*args,**kwargs)
        self.__default_vals__()
        
    def Fields(self):
        return { 'PID':None, 'Matl':None, 'T':0.0, 'NSM':0.0,
            'Z1':None, 'Z2':None }
                                                                                                                                   
class PropSolid(PropRod):
    
    def __init__(self,*args,**kwargs):
        PropRod.__init__(self,*args,**kwargs)
        self.__default_vals__()
        
    def Fields(self):
        return { 'PID':None, 'Matl':None }
                   
      