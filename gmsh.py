import numpy as np
import basic

import material as matl
import prop
import element as elem
import bcs

#import pdb

class GmshElement(object):
    
    def __init__(self, line):
        fields = line.split(' ')
        self.eid = int( fields[0] )
        self.etype = int( fields[1] )
        ntags = int( fields[2] )
        self.tags = [ ]
        self.physid = int( fields[3] )
        self.geomid = int( fields[4] )
        self.partitions = []
        i=5
        if (ntags > 2):
            npart = int( fields[i] )
            i+=1
            while i<6+npart:
                self.partitions.append( int( fields[i] ) )
                i+=1
        self.conn = []
        while i<len(fields):
            self.conn.append( int( fields[i] ) )
            i+=1
            
    def ConvertElement(self):
        if self.etype == 1:
            return elem.ElemLine2(self.conn,self.physid)
            
        elif self.etype == 2:
            return elem.ElemTria3(self.conn,self.physid)
            
        elif self.etype == 3:
            return elem.ElemQuad4(self.conn,self.physid)
            
        elif self.etype == 4:
            return elem.ElemTetra4(self.conn,self.physid)
            
        elif self.etype == 5:
            return elem.ElemElemHexa8(self.conn,self.physid)
            
        else:
            print 'Element type not supported in GmshElement.ConvertElement'
            return 0 
            
class GmshInput(object):
    
    def __init__(self,fname=None):
        self.Filename = fname
        self.PhysIDs = set()
        self.NodeSetIDs = set()
        self.SideSetIDs = set()
        
    def AddPhysicalIDs(self,ids):
        for pid in ids:
            self.PhysIDs.add(int(ids))
      
    def PhysicalIDs(self,pids):
        
        f = open( self.Filename, 'r' )
        
        e = -2
        numElem = 0
        for line in f:
            if ( line[0:9] == '$Elements' ):
                e = -1
            elif (line[0:12] == '$EndElements'):
                e = -2
            elif  e==-1:
                numElem = int(line)
                e = 0
            elif e >= 0:
                elem = GmshElement(line)
                try:
                    pids.add( elem.physid )
                except:
                    print 'Error in adding physical id to set in GmshInput.PhysicalIDs'
                    return

            if e>numElem:
                return
                    
                    
    #def FieldInt(self,i,line):
    #    try:
    #        return int( line[8*(i-1):8*i] )
    #    except:
    #        return 0
    # 
    #def FieldFloat(self,i,line):
    #    try:
    #        return float( line[8*(i-1):8*i] )
    #    except:
    #        return 0.0
   
                
    #def ReadNodes(self,nids):
    #    
    #    f = open( self.Filename, 'r' )
    #    
    #    nn = 0
    #    for line in f:
    #        if ( line[:4] == 'GRID' ):
    #            nn = nn +1
    #    
    #    self.node = np.zeros( (nn,3), basic.FLOAT_TYPE )
    #    self.nid  = np.zeros( nn, basic.INDX_TYPE )
    #    
    #    f.seek(0)
    #    n=0
    #    for line in f:
    #        line = line.ljust(80)
    #        if ( line[:4] == 'GRID' ):
    #            self.node[n,0] = self.FieldFloat(3,line)
    #            self.node[n,1] = self.FieldFloat(4,line)
    #            self.node[n,2] = self.FieldFloat(5,line)
    #            self.nid[n] = self.FieldInt(2,line)
    #            n=n+1
    #            
    #    #self.nmap = {}
    #    nn = 0
    #    for n in self.nid:
    #        self.nmap[n] = nn
    #        nn = nn + 1
        
    def ReadElements(self,elements):
        f = open( self.Filename, 'r' )
        
        if ( len(self.PhysIDs) == 0 ):
            AllPIDs = True 
        else:
            ALLPIDs = False
        
        for line in f:
            
            e = -2
            numElem = 0
            for line in f:
                if ( line[0:9] == '$Elements' ):
                    e = -1
                elif (line[0:12] == '$EndElements'):
                    e = -2
                elif  e==-1:
                    numElem = int(line)
                    e = 0
                elif e >= 0:
                    elem = GmshElement(line)
                    if ( (elem.physid in self.PhysIDs) or AllPIDs ):
                        elements.append( elem.ConvertElement() )
                    e += 1
                if e>numElem:
                    return
                    

   
    
# ---------------------
gmshfile = GmshInput('channel.msh')
pids = set()
gmshfile.PhysicalIDs(pids)  
#gmshfile.
elements=[]
gmshfile.ReadElements(elements)      
    