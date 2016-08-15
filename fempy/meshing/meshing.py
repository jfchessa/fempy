import numpy as np

#import basic
#import element as elem

#import pdb

class BiasNone(object):
    """null bias object"""
    def xi(self,n):
        return np.linspace(0,1,n)
        
class BiasExp(object):
    """exponential bias BiasExp(b)"""
    def __init__(self,b):
        self.b=b
    
    def xi(self,n):
        return (np.linspace(0,1,n))**(self.b)
        
def node_array1d(p1,p2,n,bias=BiasNone()):
    """node_array1d(p1,p2,n,bias=BiasNone())"""
    xis = bias.xi(n)
    
    pts = np.outer(np.ones(n),p1) + np.outer(xis,(p2-p1))
    
    return pts
    
        
def node_array2d(corners,n1,n2,bias1=BiasNone(),bias2=BiasNone()):
    """node_array2d(corners,n1,n2,bias1=BiasNone(),bias2=BiasNone())"""
    nn=n1*n2
    nd = corners.shape[1]
    pts = np.zeros((nn,nd))
    
    xi = 2*bias1.xi(n1)-1.0
    eta = 2*bias2.xi(n2)-1.0
    
    n=0
    for t in eta:
        for s in xi:
            N = 0.25*np.array([(1-s)*(1-t),(1+s)*(1-t), \
                    (1+s)*(1+t),(1-s)*(1+t)], dtype=basic.FLOAT_TYPE )
            pts[n,:]=np.dot(N,corners)
            n = n + 1
    
    return pts
    
        
def node_array3d(corners,n1,n2,n3,bias1=BiasNone(),bias2=BiasNone(),bias3=BiasNone()):
    """node_array3d(corners,n1,n2,n3,bias1=BiasNone(),bias2=BiasNone(),bias3=BiasNone())"""
    nn=n1*n2*n3
    nd = corners.shape[1]
    pts = np.zeros((nn,nd))
    
    xi = 2*bias1.xi(n1)-1.0
    eta = 2*bias2.xi(n2)-1.0
    zeta = 2*bias2.xi(n3)-1.0
    n=0
    eigth=float(1.0/8.0)
    for t in zeta:
        for s in eta:
            for r in xi:
                N =  eigth*np.array( [(1-r)*(1-s)*(1-t), \
                     (1+r)*(1-s)*(1-t), \
                     (1+r)*(1+s)*(1-t), \
                     (1-r)*(1+s)*(1-t), \
                     (1-r)*(1-s)*(1+t), \
                     (1+r)*(1-s)*(1+t), \
                     (1+r)*(1+s)*(1+t), \
                     (1-r)*(1+s)*(1+t) ], dtype=basic.FLOAT_TYPE )
                pts[n,:]=np.dot(N,corners)
                n = n + 1
    
    return pts
        
def gen_conn2d(ptrn,nu,nv,incu=1,incv=2):
    
    nc = ptrn.size
    conn = np.zeros( (nu*nv,nc), dtype=basic.INDX_TYPE )
    inc = np.zeros( nc, dtype=basic.INDX_TYPE )
    e=0
    for r in xrange(nv):
        for c in xrange(nu):
            conn[e] = ptrn + inc
            inc = inc + incu
            e = e + 1
        inc = inc + (incv - incu)
    
    return conn

def gen_conn3d(ptrn,nu,nv,nw,incu=1,incv=2,incw=None):
    
    if ( incw == None ):
        incw = nu+3
    
    nc = ptrn.size
    conn = np.zeros( (nu*nv*nw,nc), dtype=basic.INDX_TYPE )
    inc = np.zeros( nc, dtype=basic.INDX_TYPE )
    e=0
    for l in xrange(nw):
        for r in xrange(nv):
            for c in xrange(nu):
                conn[e] = ptrn +inc
                inc = inc +incu
                e = e + 1
            inc = inc +incv - incu
        inc = inc + incw - incv
    
    return conn
    
class MeshQuad4(object):
    
    def __init__(self,corners,n1,n2,prop,bias1=BiasNone(),bias2=BiasNone()):
        
        self.node = node_array2d(corners,n1,n2,bias1,bias2)
        self.n1 = n1
        self.n2 = n2
        
        ptrn = np.array( [0,1,n1+1,n1], dtype=basic.INDX_TYPE )
        self.conn = gen_conn2d(ptrn,n1-1,n2-1)
        
        self.element = []
        self.prop = prop
        for e in xrange( self.conn.shape[0] ):
            self.element.append( elem.ElemQuad4( self.conn[e], self.prop ) )
    
    def NumNodes(self):
        return len(self.node)
        
    def NumElem(self):
        return len(self.element)
    
    def VertexNID(self,vertid):
        """
        vertid = 0 is first corner
        vertid = 1 is second corner
        vertid = 2 is third corner
        vertid = 3 is fourth corner
        """
        if ( vertid == 0 ):
            return 0
        elif ( vertid == 1 ):
            return self.n1-1
        elif ( vertid == 2 ):
            return self.n1 * self.n2 - 1
        elif ( vertid == 3 ):
            return self.n1 * self.n2 - self.n1
            
        else:
            print "vertex id must be between 0 and 3"
            return []
            
    def EdgeNIDs(self,edgeid):
        """
        edgeid = 0 Edge between corner nodes 0 and 1
        edgeid = 1 Edge between corner nodes 1 and 2
        edgeid = 2 Edge between corner nodes 2 and 3
        edgeid = 3 Edge between corner nodes 3 and 0
        """
        if ( edgeid == 0 ):
            edgeNodes = np.array( range(0,self.n1), basic.INDX_TYPE )
        
        elif ( edgeid == 1 ):  
            edgeNodes = np.array( range(self.VertexNID(1),self.NumNodes(),self.n1), basic.INDX_TYPE )
            
        elif ( edgeid == 2 ):    
            edgeNodes = np.array( range(self.VertexNID(2),self.VertexNID(3)-1,-1), basic.INDX_TYPE )
            
        elif ( edgeid == 3 ):   
            edgeNodes = np.array( range(self.VertexNID(3),-1,-self.n1), basic.INDX_TYPE )
            
        else:
            print "Edge id must be between 0 and 3"
            return np.array( [], basic.INDX_TYPE )
            
        return edgeNodes
        
    def EdgeElem(self,edgeid):
        edgeNodes = self.EdgeNIDs(edgeid)
        edgeElem = []
        prop = self.element[0].Property()
        n0 = edgeNodes[0]
        for n1 in edgeNodes[1:]:
            edgeElem.append( elem.ElemLine2( np.array([n0,n1], basic.INDX_TYPE ), prop ))
            n0 = n1
        return edgeElem      
        
    def EdgeNormals(self,edgeid):
        edgeElem = self.EdgeElem(edgeid)
        
        if ( edgeid == 0 ):
            ninc = self.n1
        elif ( edgeid == 1 ):  
            ninc = -1
        elif ( edgeid == 2 ):    
            ninc = -self.n1
        elif ( edgeid == 3 ):   
            ninc = 1
        
        normals = []
        ve=np.zeros(3,basic.FLOAT_TYPE)
        vt=np.zeros(3,basic.FLOAT_TYPE)
        v=np.zeros(3,basic.FLOAT_TYPE)
        sdim=self.node.shape[1]
        for e in edgeElem:
            na = e.conn[0]
            nb = e.conn[1]
            nc = nb + ninc
            nd = na + ninc
            ve[:sdim] = self.node[nb]-self.node[na]
            vt[:sdim] =0.5*( self.node[nd]-self.node[na] + self.node[nc]-self.node[nb] )
            v = np.cross( ve, vt )
            enorm = np.cross(ve,v)
            enorm = enorm/np.linalg.norm(enorm)
            normals.append(enorm)
        return normals
     
    
class MeshHexa8(MeshQuad4):
    
    def __init__(self,corners,n1,n2,n3,prop,bias1=BiasNone(),bias2=BiasNone(),bias3=BiasNone()):
        
        self.node = node_array3d(corners,n1,n2,n3,bias1,bias2,bias3)
        
        dn=n1*n2
        ptrn = np.array( [0,1,n1+1,n1, dn,1+dn,n1+1+dn,n1+dn], dtype=basic.INDX_TYPE )
        self.conn = gen_conn3d(ptrn,n1-1,n2-1,n3-1)
        
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        
        self.element = []
        self.prop = prop
        for e in xrange( self.conn.shape[0] ):
            self.element.append( elem.ElemHexa8( self.conn[e], self.prop ) ) 
     
       
    def VertexNID(self,vertid):
        if ( vertid==0 ):
            return 0
        elif ( vertid==1 ):
            return self.n1-1
        elif ( vertid==2 ):
            return self.n1*self.n2-1
        elif ( vertid==3 ):
            return (self.n2-1)*self.n1
        elif (vertid==4 ):
            return (self.n3-1)*self.n1*self.n2
        elif (vertid==5 ):
            return (self.n3-1)*self.n1*self.n2+self.n1-1
        elif (vertid==6 ):
            return self.n1*self.n2*self.n3-1
        elif (vertid==7 ):
            return self.n1*self.n2*self.n3-self.n1
            
    def EdgeNIDs(self,edgeid):
        print "EdgeNIDs for MeshHexa8 not yet implemented"
        return []      
        
    def EdgeElem(self,edgeid):
        print "EdgeElem for MeshHexa8 not yet implemented"
        return []      
          
    def FaceNIDs(self,faceid):
        """faceid 0 is face in +1
        faceid 1 is face in -1
        faceid 2 is face in +2
        faceid 3 is face in -2
        faceid 4 is face in +3
        faceid 5 is face in -3"""
        if ( faceid==0 ): 
            nn = self.n2*self.n3  
            nids=np.zeros( nn, basic.INDX_TYPE )
            ii=self.VertexNID(1)
            iii=0
            for j in xrange(self.n3):
                for i in xrange(self.n2):
                    nids[iii] = ii
                    ii=ii+self.n1
                    iii=iii+1
                #ii = ii + 0
                
        elif ( faceid==1 ):
            nn = self.n2*self.n3
            nids=np.zeros( nn, basic.INDX_TYPE )
            nn = self.n2*self.n3  
            nids=np.zeros( nn, basic.INDX_TYPE )
            ii=self.VertexNID(0)
            iii=0
            for j in xrange(self.n3):
                for i in xrange(self.n2):
                    nids[iii] = ii
                    ii=ii+self.n1
                    iii=iii+1
                #ii = ii + 0
                
        elif ( faceid==2 ):
            nn = self.n1*self.n3
            nids=np.zeros( nn, basic.INDX_TYPE )
            ii=self.VertexNID(3)
            iii=0
            for j in xrange(self.n3):
                for i in xrange(self.n1):
                    nids[iii] = ii
                    ii=ii+1
                    iii=iii+1
                ii = ii + self.n1*(self.n2-1)
                
        elif ( faceid==3 ):
            nn = self.n1*self.n3
            nids=np.zeros( nn, basic.INDX_TYPE )
            ii=self.VertexNID(0)
            iii=0
            for j in xrange(self.n3):
                for i in xrange(self.n1):
                    nids[iii] = ii
                    ii=ii+1
                    iii=iii+1
                ii = ii + self.n1*(self.n2-1)
                
        elif ( faceid==4 ):
            nn = self.n1*self.n2
            nids=np.zeros( nn, basic.INDX_TYPE )
            ii=self.VertexNID(4)
            iii=0
            for j in xrange(self.n2):
                for i in xrange(self.n1):
                    nids[iii] = ii
                    ii=ii+1
                    iii=iii+1
                
        elif ( faceid==5 ):
            nn = self.n1*self.n2
            nids=np.zeros( nn, basic.INDX_TYPE )
            ii=self.VertexNID(0)
            iii=0
            for j in xrange(self.n2):
                for i in xrange(self.n1):
                    nids[iii] = ii
                    ii=ii+1
                    iii=iii+1
            
        return nids
        
    def FaceElem(self,faceid):
        
        if ( faceid==0 ):
            n0 = self.VertexNID(1)
            dz=self.n1*self.n2
            ptrn = np.array( [n0,n0+self.n1,n0+dz+self.n1,n0+dz], dtype=basic.INDX_TYPE )
            fconn = gen_conn2d(ptrn,self.n2-1,self.n3-1,self.n1,2*(self.n1))
           
        elif ( faceid==1 ):
            n0 = self.VertexNID(0)
            dz=self.n1*self.n2
            ptrn = np.array( [n0,n0+self.n1,n0+dz+self.n1,n0+dz], dtype=basic.INDX_TYPE )
            fconn = gen_conn2d(ptrn,self.n2-1,self.n3-1,self.n1,2*(self.n1))
           
        elif ( faceid==2 ):
            n0 = self.VertexNID(3)
            dz=self.n1*self.n3
            ptrn = np.array( [n0,n0+1,n0+dz+1,n0+dz], dtype=basic.INDX_TYPE )
            fconn = gen_conn2d(ptrn,self.n1-1,self.n3-1,1,dz-self.n1+2)
           
        elif ( faceid==3 ):
            n0 = self.VertexNID(0)
            dz=self.n1*self.n3
            ptrn = np.array( [n0,n0+1,n0+dz+1,n0+dz], dtype=basic.INDX_TYPE )
            fconn = gen_conn2d(ptrn,self.n1-1,self.n3-1,1,dz-self.n1+2)
           
        elif ( faceid==4 ):
            n0 = self.VertexNID(4)
            ptrn = np.array( [n0,n0+1,n0+self.n1+1,n0+self.n1], dtype=basic.INDX_TYPE )
            fconn = gen_conn2d(ptrn,self.n1-1,self.n2-1,1,2)
           
        elif ( faceid==5 ):
            n0 = self.VertexNID(0)
            ptrn = np.array( [n0,n0+1,n0+self.n1+1,n0+self.n1], dtype=basic.INDX_TYPE )
            fconn = gen_conn2d(ptrn,self.n1-1,self.n2-1,1,2)
        
        felem = []
        for ec in fconn:
            felem.append( elem.ElemQuad4( ec, self.prop ) )
            
        return felem
    
    def FaceNormals(self,faceid):
        
        felem = self.FaceElem(faceid)
        fnorm=[]
        for f in felem:
            (na,nb,nc,nd) = f.conn
            v1 = self.node[nc]-self.node[na]
            v2 = self.node[nd]-self.node[nb]
            fn = np.cross(v1,v2)
            fn = fn/np.linalg.norm(fn)
            fnorm.append(fn)
        return fnorm
                               