import numpy as np

import basic

def quadrature_gauss1d(numpt):
    
    if ( numpt==1 ):
        quadpoint = np.array( [0.000000000000000], dtype = basic.FLOAT_TYPE )
        quadweight = np.array( [2.000000000000000], dtype = basic.FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 2 ):
        quadpoint = np.array( [0.577350269189626, -0.577350269189626], dtype = basic.FLOAT_TYPE )
        quadweight = np.array( [1.0, 1.0 ], dtype = basic.FLOAT_TYPE )  
        return [quadpoint,quadweight]
        
    if ( numpt == 3 ):
        quadpoint = np.array( [0.774596669241483, -0.774596669241483, 0.000000000000000], dtype = basic.FLOAT_TYPE )
        
        quadweight = np.array( [0.555555555555556, 0.555555555555556, 0.888888888888889], dtype = basic.FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 4 ):
        quadpoint = np.array( [0.861134311594053,-0.861134311594053, \
                0.339981043584856,-0.339981043584856], dtype = basic.FLOAT_TYPE )
        
        quadweight = np.array( [0.347854845137454, 0.347854845137454,\
                0.652145154862546, 0.652145154862546], dtype = basic.FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 5 ):
        quadpoint = np.array( [0.906179845938664,-0.906179845938664,0.538469310105683,\
             -0.538469310105683, 0.000000000000000], dtype = basic.FLOAT_TYPE )
        
        quadweight = np.array( [0.236926885056189, 0.236926885056189, 0.478628670499366,\
                0.478628670499366, 0.568888888888889], dtype = basic.FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 6 ):
        quadpoint = np.array( [0.932469514203152,-0.932469514203152,\
            0.661209386466265,-0.661209386466265,0.238619186003152,\
            -0.238619186003152], dtype = basic.FLOAT_TYPE )
        
        quadweight = np.array( [0.171324492379170,0.171324492379170,0.360761573048139,\
                0.360761573048139,0.467913934572691, 0.467913934572691], dtype = basic.FLOAT_TYPE )   
        return [quadpoint,quadweight]
        
    if ( numpt == 7 ):
        quadpoint = np.array( [0.949107912342759,-0.949107912342759,0.741531185599394,\
                            -0.741531185599394,0.405845151377397,\
                            -0.405845151377397,0.000000000000000], dtype = basic.FLOAT_TYPE )
    
        quadweight = np.array( [0.129484966168870,0.129484966168870,0.279705391489277,\
                0.279705391489277,0.381830050505119,0.381830050505119,\
                0.417959183673469 ], dtype = basic.FLOAT_TYPE )  
        return [quadpoint,quadweight]
        
    if ( numpt == 8 ):
        quadpoint = np.array( [0.960289856497536,-0.960289856497536,0.796666477413627,\
         -0.796666477413627,0.525532409916329,-0.525532409916329,0.183434642495650,\
         -0.183434642495650], dtype = basic.FLOAT_TYPE )
        
        quadweight = np.array( [0.101228536290376, 0.101228536290376, 0.222381034453374,\
            0.222381034453374, 0.313706645877887, 0.313706645877887,\
            0.362683783378362, 0.362683783378362 ], dtype = basic.FLOAT_TYPE )  
        return [quadpoint,quadweight]

#def quadrature_gauss( sdim, p ):
#    
#    nn = int(np.ceil(0.5*(p+1)))
#    numpt = nn ** sdim
#    if ( sdim == 1 ):
#        return quadrature_gauss1d( numpt )
#
#    qrule = quadrature_gauss1d(nn)
#
#    quadweight = np.zeros( numpt, dtype = basic.FLOAT_TYPE )
#    quadpoint = np.zeros( (numpt,sdim), dtype = basic.FLOAT_TYPE ) 
#    
#    for s in xrange(sdim):
#        quadweight[s*nn:(s+1)*nn]  = qrule[1]
#        quadpoint[s*nn:(s+1)*nn,0] = qrule[0]
#        quadpoint[s:numpt:nn,1] = qrule[0]
#        
#    return [quadpoint,quadweight]
            
def compound_quadrature(p1,w1,p2,w2,p3=None,w3=None):
    
    if ( p3==None ): 
        
        npts = w1.size * w2.size
        sdim = 2 
        qwts = np.zeros( npts, dtype = basic.FLOAT_TYPE )
        qpts = np.zeros((npts,sdim), dtype = basic.FLOAT_TYPE )
        
        n=0
        for j in xrange(w2.size):        
            for i in xrange(w1.size):  
                qpts[n,:]  =  [ p1[i], p2[j] ]  
                qwts[n] = w1[i]*w2[j]
                n = n+1
        
    else:
        
        npts = w1.size * w2.size * w3.size
        sdim = 3
        qwts = np.zeros( npts, dtype = basic.FLOAT_TYPE )
        qpts = np.zeros((npts,sdim), dtype = basic.FLOAT_TYPE )
        
        n=0
        for k in xrange(w3.size):
            for j in xrange(w2.size):        
                for i in xrange(w1.size):  
                    qpts[n,:]  =  [ p1[i], p2[j], p3[k] ]  
                    qwts[n] = w1[i]*w2[j]*w3[k]
                    n = n+1
                          
    return [qpts,qwts]
                                                                      
def quadrature_simplex(sdim,quadorder):
    sixth = 1.0/6
    
    if ( sdim == 3 ):  # tetrahedra
#      if ( quadorder ~= 1 &  quadorder ~= 2 &  quadorder ~= 3  ) 
#        % check for valid quadrature order
#        disp('Incorect quadrature order for QUADRATURE_SIMPLEX');
#        quadorder = 1;
#       end
        
        if  ( quadorder == 1 ):
            quadpoint = np.array([[ 0.25, 0.25, 0.25 ]],dtype=basic.FLOAT_TYPE)
            quadweight = np.array([ 1.0 ],dtype=basic.FLOAT_TYPE)
            return [quadpoint,sixth*quadweight]
        
        if ( quadorder == 2 ): 
            quadpoint = np.array([[ 0.58541020,  0.13819660,  0.13819660],\
                      [ 0.13819660,  0.58541020,  0.13819660],\
                      [ 0.13819660,  0.13819660,  0.58541020],\
                      [ 0.13819660,  0.13819660,  0.13819660]],dtype=basic.FLOAT_TYPE)
            quadweight = np.array([.25,.25,.25,.25],dtype=basic.FLOAT_TYPE)
            return [quadpoint,sixth*quadweight]
        
        if ( quadorder == 3 ):
            quadpoint = np.array( [[ 0.25,  0.25,  0.25],\
                     [ 0.5,   sixth,   sixth ],\
                     [ sixth,   0.5,   sixth ],\
                     [ sixth,   sixth,   0.5 ],\
                     [ sixth,   sixth,   sixth ]], dtype=basic.FLOAT_TYPE )
            quadweight = np.array( [-0.8, .45, .45, .45, .45], dtype=basic.FLOAT_TYPE )
            return [quadpoint,sixth*quadweight]
        
         
    else:  # TRIANGLES
      
#      if ( quadorder > 7 ) % check for valid quadrature order
#        disp('Quadrature order too high for QUADRATURE_SIMPLEX');
#        quadorder = 1;
#      end
      
        if ( quadorder <= 1 ):  # set quad points and quadweights
            quadpoint = np.array([[ 0.3333333333333, 0.3333333333333 ]], dtype=basic.FLOAT_TYPE )
            quadweight =np.array([ 0.5 ], dtype=basic.FLOAT_TYPE )
            return [quadpoint,quadweight]  
        
        if ( quadorder == 2 ): 
            quadweight = sixth*np.ones( 3, dtype=basic.FLOAT_TYPE )
            quadpoint = np.array( [[ 0.1666666666667, 0.1666666666667 ],\
                            [ 0.6666666666667, 0.1666666666667 ],\
                            [ 0.1666666666667, 0.6666666666667 ]], dtype=basic.FLOAT_TYPE ) 
            return [quadpoint,quadweight]  
        
        if ( quadorder <= 5 ): 
        
            quadpoint = np.array( [[ 0.1012865073235, 0.1012865073235 ],\
                        [ 0.7974269853531, 0.1012865073235 ],\
                        [ 0.1012865073235, 0.7974269853531 ],\
                        [ 0.4701420641051, 0.0597158717898 ],\
                        [ 0.4701420641051, 0.4701420641051 ],\
                        [ 0.0597158717898, 0.4701420641051 ],\
                        [ 0.3333333333333, 0.3333333333333 ]], dtype=basic.FLOAT_TYPE )
        
            quadweight = np.array( [0.1259391805448, 0.1259391805448, 0.1259391805448,\
                    0.1323941527888, 0.1323941527885,  0.1323941527885, \
                    0.2250000000000 ], dtype=basic.FLOAT_TYPE )
                    
            return [quadpoint,0.5*quadweight]    
        
        quadpoint = np.array( [[ 0.0651301029022, 0.0651301029022 ],\
                    [ 0.8697397941956, 0.0651301029022 ],\
                    [ 0.0651301029022, 0.8697397941956 ],\
                    [ 0.3128654960049, 0.0486903154253 ],\
                    [ 0.6384441885698, 0.3128654960049 ],\
                    [ 0.0486903154253, 0.6384441885698 ],\
                    [ 0.6384441885698, 0.0486903154253 ],\
                    [ 0.3128654960049, 0.6384441885698 ],\
                    [ 0.0486903154253, 0.3128654960049 ],\
                    [ 0.2603459660790, 0.2603459660790 ],\
                    [ 0.4793080678419, 0.2603459660790 ],\
                    [ 0.2603459660790, 0.4793080678419 ],\
                    [ 0.3333333333333, 0.3333333333333 ]], dtype=basic.FLOAT_TYPE )
        
        quadweight = np.array( [ 0.0533472356088, 0.0533472356088, 0.0533472356088,\
                    0.0771137608903, 0.0771137608903, 0.0771137608903, 0.0771137608903,\
                    0.0771137608903, 0.0771137608903, 0.1756152576332, 0.1756152576332, \
                    0.1756152576332, -0.1495700444677 ], dtype=basic.FLOAT_TYPE )
                     
        return [quadpoint,0.5*quadweight]  
        
    