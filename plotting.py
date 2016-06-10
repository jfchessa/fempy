import numpy as np
import basic
import mayavi.mlab as mlab 
#import tvtk.api as tvtk
from tvtk.api import tvtk

import pdb

def make_vtk_data(node,element,ndata=None,cdata=None):
    
    numElem = len(element)
    numNode = len(node)
    
    # generate the VTK cell data.  This is a linear array of integers
    # where the first number in a pattern is the number of nodes in the 
    # connectivity followed by the connectivity
    cellSize = 0  # find the size of cells
    for e in element:
        cellSize = cellSize + e.NumNodes() + 1
        
    cell = np.zeros( cellSize, dtype=basic.INDX_TYPE )
    offset = np.zeros( numElem, dtype=basic.INDX_TYPE )
    cell_type = np.zeros( numElem, dtype=basic.INDX_TYPE )
    
    n = 1
    en = 0
    for e in element:
       
        nn = e.NumNodes()
        cell[n-1] = nn
        cell[n:n+nn] = e.conn
       
        offset[en]= n-1
                
        if ( e.Type() == 'LINE2' ):
            cell_type[en] = 3
        elif ( e.Type() == 'LINE3' ): 
            cell_type[en] = 21  
        elif ( e.Type() == 'TRIA3' ): 
            cell_type[en] = 5
        elif ( e.Type() == 'TRIA6' ): 
            cell_type[en] = 22 
        elif ( e.Type() == 'QUAD4' ): 
            cell_type[en] = 9
        elif ( e.Type() == 'QUAD8' ): 
            cell_type[en] = 23
        elif ( e.Type() == 'TETRA4' ): 
            cell_type[en] = 10
        elif ( e.Type() == 'TETRA10' ): 
            cell_type[en] = 24
        elif ( e.Type() == 'HEXA8' ): 
            cell_type[en] = 12
        elif ( e.Type() == 'HEXA20' ): 
            cell_type[en] = 25
            
        n = n + nn + 1
        en = en + 1
       
    cell_array = tvtk.CellArray()
    cell_array.set_cells(numElem, cell)
    
    points = np.zeros( (numNode,3), dtype=basic.FLOAT_TYPE )
    points[:,:node.shape[1]] = node
    
    ug = tvtk.UnstructuredGrid(points=points)
    
    # Now just set the cell types and reuse the ug locations and cells.
    ug.set_cells(cell_type, offset, cell_array)
    
    scalars = node[:,0]
    ug.point_data.scalars = scalars
    ug.point_data.scalars.name = 'scalars'

    return ug     
        
# --------------------------------------------------------------------------- #
# set up the fe mesh and the problem parameters
import prop
import material as matl
import meshing as mesh

plateWidth = 10.0
plateHeight = 10.0
he = 0.1

mat0 = matl.LinearElasticMat( 10e6, .33, 2.45e-4)
prop0 = prop.PlaneStress( mat0, thk = 1.0 )

nnx = int( plateWidth/he ) + 1
nny = int( plateHeight/he ) + 1
grid = mesh.MeshQuad4( np.array([[0,0],[plateWidth,0],[plateWidth,plateHeight],\
                 [0,plateHeight]], basic.FLOAT_TYPE ), nnx, nny, prop0 )

ugdata = make_vtk_data(grid.node,grid.element)  

from mayavi import mlab
mlab.clf()
mlab.pipeline.surface(ugdata)