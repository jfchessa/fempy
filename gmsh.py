import numpy as np
import basic

import material as matl
import prop
import element as elem
import bcs

import pdb

class GmshInput(object):
    
    def __init__(self,fname=None):
        self.filename = fname
        self.renumbered = False
        
        self.node = np.array([],basic.FLOAT_TYPE)
        self.nmap={}
        self.element = []
        self.materials = {}
        self.properties = {}
        self.forces = []
        self.spcs = []
     
    def FieldInt(self,i,line):
        try:
            return int( line[8*(i-1):8*i] )
        except:
            return 0
     
    def FieldFloat(self,i,line):
        try:
            return float( line[8*(i-1):8*i] )
        except:
            return 0.0
   
                
    def ReadNodes(self):
        
        f = open( self.filename, 'r' )
        
        nn = 0
        for line in f:
            if ( line[:4] == 'GRID' ):
                nn = nn +1
        
        self.node = np.zeros( (nn,3), basic.FLOAT_TYPE )
        self.nid  = np.zeros( nn, basic.INDX_TYPE )
        
        f.seek(0)
        n=0
        for line in f:
            line = line.ljust(80)
            if ( line[:4] == 'GRID' ):
                self.node[n,0] = self.FieldFloat(3,line)
                self.node[n,1] = self.FieldFloat(4,line)
                self.node[n,2] = self.FieldFloat(5,line)
                self.nid[n] = self.FieldInt(2,line)
                n=n+1
                
        #self.nmap = {}
        nn = 0
        for n in self.nid:
            self.nmap[n] = nn
            nn = nn + 1
        
    def ReadElements(self):
        f = open( self.filename, 'r' )
        
        #self.element=[]
        for line in f:
            line = line.ljust(80)
            if ( line[:4] == 'CROD' ):
                eid = self.FieldInt(2,line)
                pid = self.FieldInt(3,line)
                n1  = self.FieldInt(4,line)
                n2  = self.FieldInt(5,line)
                self.element.append( elem.ElemLine2( \
                    np.array( [n1,n2],basic.INDX_TYPE), prop=pid ) )
                      
            elif ( line[:4] == 'CBAR' ):
                eid = self.FieldInt(2,line)
                pid = self.FieldInt(3,line)
                n1  = self.FieldInt(4,line)
                n2  = self.FieldInt(5,line)
                self.element.append( elem.ElemLine2( \
                    np.array( [n1,n2],basic.INDX_TYPE), prop=pid ) )
                v1 = self.FieldFloat(6,line)
                v2 = self.FieldFloat(7,line)
                v3 = self.FieldFloat(8,line)
                      
            elif ( line[:6] == 'CTRIA3' ):
                eid = self.FieldInt(2,line)
                pid = self.FieldInt(3,line)
                n1  = self.FieldInt(4,line)
                n2  = self.FieldInt(5,line)
                n3  = self.FieldInt(6,line)
                self.element.append( elem.ElemTria3( \
                    np.array([n1,n2,n3],basic.INDX_TYPE), prop=pid ) )
                      
            elif ( line[:6] == 'CQUAD4' ):
                eid = self.FieldInt(2,line)
                pid = self.FieldInt(3,line)
                n1  = self.FieldInt(4,line)
                n2  = self.FieldInt(5,line)
                n3  = self.FieldInt(6,line)
                n4  = self.FieldInt(7,line)
                self.element.append( elem.ElemQuad4( np.array( \
                      [n1,n2,n3,n4],basic.INDX_TYPE), prop=pid ) )
                            
            elif ( line[:6] == 'CTETRA' ):
                eid = self.FieldInt(2,line)
                pid = self.FieldInt(3,line)
                n1  = self.FieldInt(4,line)
                n2  = self.FieldInt(5,line)
                n3  = self.FieldInt(6,line)
                n4  = self.FieldInt(7,line)
                self.element.append( elem.ElemTetra4(np.array( \
                      [n1,n2,n3,n4],basic.INDX_TYPE), prop=pid ) )
                      
            elif ( line[:5] == 'CHEXA' ):
                line = line + f.next()
                l2 = f.next().ljust(80)[8:]
                line = line+l2
                eid = self.FieldInt(2,line)
                pid = self.FieldInt(3,line)
                n1  = self.FieldInt(4,line)
                n2  = self.FieldInt(5,line)
                n3  = self.FieldInt(6,line)
                n4  = self.FieldInt(7,line)
                n5  = self.FieldInt(8,line)
                n6  = self.FieldInt(9,line)
                n7  = self.FieldInt(10,line)
                n8  = self.FieldInt(11,line)
                self.element.append( elem.ElemHexa8( np.array( \
                      [n1,n2,n3,n4,n5,n6,n7,n8],basic.INDX_TYPE), \
                       prop=pid ) )      


   
    
# ---------------------
#nasfile = NastranInput('quad4.bdf')
#nasfile.ReadMaterials()
#nasfile.ReadProperties()
#nasfile.ReadNodes()
#nasfile.ReadElements()
#nasfile.ReadSPCs()
#nasfile.ReadForces()
#
#nasfile.Renumber()
#        
    