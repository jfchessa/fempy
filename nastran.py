import numpy as np
import basic

import material as matl
import prop
import element as elem
import bcs

import pdb

class NastranInput(object):
    
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
     
    def ReadMaterials(self):
        #self.materials = {}
        
        f = open( self.filename, 'r' )
        for line in f:
            line = line.ljust(80)
            if ( line[:4] == 'MAT1' ):
                mid = self.FieldInt(2,line)
                E = self.FieldFloat(3,line)
                G = self.FieldFloat(4,line)
                nu = self.FieldFloat(5,line)
                rho = self.FieldFloat(6,line)
                if ( nu==0 ):
                    nu = 1.0 - E/(2*G)
                self.materials[mid] =  matl.LinearElasticMat( E,nu,rho )
                
    def ReadProperties(self):
        #self.properties = {}
        
        f = open( self.filename, 'r' )
        for line in f:
            line = line.ljust(80)
            
            if ( line[:4] == 'PROD' ):
                pid =  self.FieldInt(2,line)
                mid =  self.FieldInt(3,line)
                A = self.FieldFloat(4,line)
                J = self.FieldFloat(5,line)
                c = self.FieldFloat(6,line)
                self.properties[pid] = prop.Rod(mid,A,J,c)
                
            elif ( line[:4] == 'PBAR' ):
                pid =  self.FieldInt(2,line)
                mid =  self.FieldInt(3,line)
                A =  self.FieldFloat(4,line)
                I1 = self.FieldFloat(5,line)
                I2 = self.FieldFloat(6,line)
                J =  self.FieldFloat(7,line)
                self.properties[pid] = prop.Bar(mid,A,I1,I2,J)
                
            elif ( line[:6] == 'PSHELL' ):
                pid =  self.FieldInt(2,line)
                mid =  self.FieldInt(3,line)
                t =    self.FieldFloat(4,line)
                mid2 = self.FieldFloat(5,line)
                Ifac = self.FieldFloat(6,line)
                mid3 = self.FieldFloat(7,line)
                sfac = self.FieldFloat(8,line)
                self.properties[pid] = prop.Shell(mid,t,Ifac,sfac)
                
            elif ( line[:4] == 'PSOLID' ):
                pid =  self.FieldInt(2,line)
                mid =  self.FieldInt(3,line)
                
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

    def ReadForces(self):
        f = open( self.filename, 'r' )
        
        #self.forces=[]
        for line in f:
            line = line.ljust(80)
            if ( line[:5] == 'FORCE' ):
                sid = self.FieldInt(2,line)
                nid = self.FieldInt(3,line)
                cid = self.FieldInt(4,line)
                fmag = self.FieldFloat(5,line)
                n1 = self.FieldFloat(6,line)
                n2 = self.FieldFloat(7,line)
                n3 = self.FieldFloat(8,line)
                fvec = fmag*np.array([n1,n2,n3],basic.FLOAT_TYPE)
                self.forces.append( bcs.Force(nid,fvec) )
                
    def CID2LDOF(self,cid):
        ldof=[]
        cid = cid.strip()
        if ( len(cid) == 0 ):
            return ldof
        for c in cid:
            ldof.append( int(c)-1 )
        return ldof

    def ReadSPCs(self):
        f = open( self.filename, 'r' )
        
        #self.spcs = []
        for line in f:
            line = line.ljust(80)
            if ( line[:3] == 'SPC' ):
                sid = self.FieldInt(2,line)
                nid1 = self.FieldInt(3,line)
                c1  = line[24:32]
                d1  = self.FieldFloat(5,line)
                nid2 = self.FieldInt(6,line)
                c2  = line[56:64]
                d2  = self.FieldFloat(8,line)
                
                for s in self.CID2LDOF(c1):
                    self.spcs.append( bcs.SPC( nid1,s,d1) )
                for s in self.CID2LDOF(c2):
                    self.spcs.append( bcs.SPC( nid2,s,d2) )         
        
    def Renumber(self):
        
        if ( self.renumbered ):
            return
            
        #pdb.set_trace()
        # set materials in properties
        for pid, prop in self.properties.iteritems():
            mid = int(prop.material)
            try:
                prop.material = self.materials[mid]
            except:
                print 'Invalid material id in Renumber()'
               
        # set properties in elements
        for e in self.element:
            pid = int(e.prop)
            try:
                e.prop = self.properties[pid]
            except:
                print 'invalid property id in Renumber()'
        
        # renumber element connectivity
        for e in self.element:
            #pdb.set_trace()
            n=0
            for i in e.conn:
                try:
                    e.conn[n] = self.nmap[i]
                except:
                    print 'Invalid node id in Renumber()'
                n=n+1
        
        # renumber force
        for f in self.forces:
            try:
                f.nid = self.nmap[ f.nid ]
            except:
                print 'Invalid node id in Renumber()'
        
        # renumber spc
        for s in self.spcs:
            try:
                s.nid = self.nmap[ s.nid ]
            except:
                print 'Invalid node id in Renumber()'
        
        self.renumbered = True                
                                            
# ---------------------
nasfile = NastranInput('quad4.bdf')
nasfile.ReadMaterials()
nasfile.ReadProperties()
nasfile.ReadNodes()
nasfile.ReadElements()
nasfile.ReadSPCs()
nasfile.ReadForces()

nasfile.Renumber()
        
    