import numpy as np
#import basic
#import element as elem
import fempy as fp

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
        self.nidmap=dict()
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
            return elem.ElemHexa8(self.conn,self.physid)
            
        else:
            print 'Element type not supported in GmshElement.ConvertElement'
            return 0 
            
class GmshInput(object):
    
    def __init__(self,fname=None):
        self.Filename = fname
        
        self.PhysIDs = set()
        self.NodeSetIDs = set()
        self.SideSetIDs = set()
        
#        self.Renumber = True
#        self.nidmap = None
        
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
   
        #if self.nidmap == None:
        #    self.ReadNodeMap()
            
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
                        #elem.conn = self.RenumberNIDs(elem.conn)
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
                        #elem.conn = self.RenumberNIDs(elem.conn)
                        sidesets[elem.physid].append(elem.ConvertElement())
                    e += 1
                if e>numElem:
                    break
            
            return sidesets
                                                      
    #def ReadNodeMap(self):
    #         
    #    f = open( self.Filename, 'r' )
    #  
    #    for line in f:
    #        
    #        n = -2
    #        numNode = 0
    #        for line in f:
    #            
    #            if ( line[0:6] == '$Nodes' ):
    #                n = -1
    #                
    #            elif (line[0:9] == '$EndNodes'):
    #                n = -2
    #                
    #            elif  n==-1:
    #                numNode = int(line)
    #                nids = {}
    #                
    #                if not self.Renumber:
    #                    for i in xrange(numNode):
    #                        nids[i]=i
    #                    return nids
    #                    
    #                n = 0
    #    
    #            elif n >= 0:
    #                fields = line.split(' ')
    #                i = int(fields[0])
    #                nids[i] = n
    #                n += 1
    #            if n>numNode:
    #                break
    #          
    #    if self.Renumber:
    #        self.nidmap = nids   
    #                   
    #    return nids    
                                                              
    def ReadNodes(self):
       
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
                    nodes = basic.NodeArray() #np.zeros((numNode,3),float)
                    #nids = {}
                    n = 0
                elif n >= 0:
                    fields = line.split(' ')
                    i = int(fields[0])
                    x = float(fields[1])
                    y = float(fields[2])
                    z = float(fields[3])
                    #nids[i] = n
                    nodes[i] = np.array([x, y, z])
                    n += 1
                if n>numNode:
                    break
             
            #if self.Renumber:
            #    self.nidmap = nids
                               
            #return [nodes,nids]
            return nodes
                    
    def ReadElements(self,propMap=None):
            
        f = open( self.Filename, 'r' )
        elements = elem.ElementArray() #[]
        
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
                    elm = GmshElement(line)
                    if ( (elm.physid in self.PhysIDs) or AllPIDs ):
                        #self.RenumberNIDs(elem.conn)
                        elements[elm.eid] =  elm.ConvertElement() 
                    e += 1
                if e>numElem:
                    break
            
        if ( not propMap == None ):
            for eid, e in elements.iteritems():
                try:
                    e.prop = propMap[e.prop]
                except:
                    a=1
                        
        return elements
     
    #def RenumberNIDs(self,nids):
    #    if not self.Renumber:
    #        return
    #        
    #    if self.nidmap == None:
    #        self.ReadNodeMap()
    #        
    #    for i in xrange(len(nids) ):
    #        ni=nids[i]
    #        try:
    #            nids[i] = self.nidmap[ni]  
    #        except KeyError:
    #            print 'nid ', i, ' not found in node map for RenumberNIDs()'
    #    return nids
            
    #def RenumberConnectivity(self,nidmap,elements):
    #    
    #    for e in elements:
    #        e.conn = self.RenumberNIDs(e.conn)
    #    
    #    return elements
    #    
    #def RenumberNodeSets(self,nidmap,nsets):
    #    
    #    if self.nidmap == None:
    #        self.ReadNodeMap()
    #        
    #    for pid, ns in nsets.iteritems():
    #        news = set()
    #        while ns:
    #            i = ns.pop()
    #            news.add( nidmap[i] )
    #        nsets[pid]=news

#-----------------------------------------------------------------------------
#gmshfile = GmshInput('channel.msh')
#
#channelPID = 163
#fluidPID = 162
#
#inletPID = 164
#heatFluxPID = 165
#
##gmshfile.AddPhysicalIDs([channelPID,fluidPID])
#gmshfile.AddPhysicalIDs([channelPID])
#
#gmshfile.AddSideSetIDs( inletPID )
#gmshfile.AddNodeSetIDs( heatFluxPID )
# 
#nodes = gmshfile.ReadNodes()  
#elements = gmshfile.ReadElements()  
#    
#nodesets = gmshfile.ReadNodeSets()
#sidesets = gmshfile.ReadSideSets()
#           