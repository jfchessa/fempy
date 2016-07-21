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
        try:
            for pid in ids:
                self.PhysIDs.add(int(pid))
        except TypeError:
            self.PhysIDs.add(int(ids))
      
    def AddSideSetIDs(self,ids):
        try:
            for pid in ids:
                self.SideSetIDs.add(int(pid))
        except TypeError:
            self.SideSetIDs.add(int(ids))
            
    def AddNodeSetIDs(self,ids):
        try:
            for pid in ids:
                self.NodeSetIDs.add(int(pid))
        except TypeError:
            self.NodeSetIDs.add(int(ids))
            
      
    def PhysicalIDs(self):
        
        pids = set()
        
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

            if e>numElem:
                break
        return pids
    
    def NumNodes(self):
       
        f = open( self.Filename, 'r' )
      
        for line in f:
            
            n = -2
            for line in f:
                if ( line[0:6] == '$Nodes' ):
                    n = -1
                elif (line[0:9] == '$EndNodes'):
                    n = -2
                elif  n==-1:
                    return int(line)

    def ReadNodeSets(self):
   
        nodesets = {}
        for pid in self.NodeSetIDs:
            nodesets[pid] = set()
        
        f = open( self.Filename, 'r' )
        
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
                    if ( (elem.physid in self.NodeSetIDs) ):
                        for i in elem.conn:
                            nodesets[elem.physid].add( i )
                    e += 1
                if e>numElem:
                    break
            return nodesets
                                               
    def ReadSideSets(self):
   
        sidesets = {}
        for pid in self.SideSetIDs:
            sidesets[pid] = []
        
        f = open( self.Filename, 'r' )
        
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
                    if ( (elem.physid in self.SideSetIDs) ):
                        sidesets[elem.physid].append(elem.ConvertElement())
                    e += 1
                if e>numElem:
                    break
            
            return sidesets
                                                       
    def ReadNodes(self,nodes,nids):
       
        f = open( self.Filename, 'r' )
      
        for line in f:
            
            n = -2
            numNode = 0
            for line in f:
                if ( line[0:6] == '$Nodes' ):
                    n = -1
                elif (line[0:9] == '$EndNodes'):
                    n = -2
                elif  n==-1:
                    numNode = int(line)
                    n = 0
                elif n >= 0:
                    fields = line.split(' ')
                    i = int(fields[0])
                    x = float(fields[1])
                    y = float(fields[2])
                    z = float(fields[3])
                    nids[i] = n
                    nodes[n] = np.array([x, y, z])
                    n += 1
                if n>numNode:
                    break
                          
        
    def ReadElements(self,elements):
        f = open( self.Filename, 'r' )
        
        if ( len(self.PhysIDs) == 0 ):
            AllPIDs = True 
        else:
            AllPIDs = False
        
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
                    break
                    
# -----------------------------------
gmshfile = GmshInput('channel.msh')
pids = gmshfile.PhysicalIDs()  
gmshfile.AddPhysicalIDs([16,17])

elements=[]
gmshfile.ReadElements(elements)  

nn = gmshfile.NumNodes()
nodes=np.zeros((nn,3),float)
nids = { }
gmshfile.ReadNodes(nodes,nids)    
    
gmshfile.AddSideSetIDs( [19,20] )
gmshfile.AddNodeSetIDs( 18 )
nodesets = gmshfile.ReadNodeSets()
sidesets = gmshfile.ReadSideSets()
    