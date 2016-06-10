import numpy as np

class Rod(object):
    
    def __init__(self,mat,A,J=0.0,c=0.0):
        self.material = mat
        self.area = A
        self.torsionStiff = J
        self.torsionRec = c
        
    def __repr__(self):
        return 'Rod property card:\n'+'A='+str(self.area)+' J='+\
            str(self.torsionStiff)+' c='+str(self.torsionRec)+\
            str(self.material)

class Bar(Rod):
    
    def __init__(self,mat,A,I1,I2,J):
        self.material = mat
        self.area = A
        self.torsionStiff = J
        self.I1 = I1
        self.I2 = I2
        
    def __repr__(self):
        return 'Bar property card:\n'+'A='
    
        
class PlaneStress(Rod):
    
    def __init__(self,mat,thk=1.0):
        self.material = mat
        self.thk = float(thk)
        
    def __repr__(self):
        return 'Plane Stress property card: t='+str(self.thk)+'\n'+str(self.material)       
    
class Shell(PlaneStress):
    
    def __init__(self,mat,thk=1.0,ifac=1.0,sfac=1.0):
        self.material = mat
        self.thk = float(thk)
        self.ifac = float(ifac)
        self.sfac = float(sfac)
        
    def __repr__(self):
        return 'Shell property card: t='+str(self.thk)+'\n'+str(self.material)       
    
class PlaneStrain(PlaneStress):
    
    def __init__(self,mat):
        self.material = mat
        
    def __repr__(self):
        return 'Shell property card:\n'+str(self.material)
        
class Solid3D(PlaneStress):
    
    def __init__(self,mat):
        self.material = mat
        
    def __repr__(self):
        return '3D solid property card:\n'+str(self.material)
                         
class Axisymmetric(PlaneStress):
      
    def __init__(self,mat):
        self.material = mat
        
    def __repr__(self):
        return 'Axisymmetric property card:\n'+str(self.material)
