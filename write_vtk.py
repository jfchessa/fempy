import numpy as np
import basic

import pdb
# pdb.set_trace()

class OutputFile(object):
    
    def __init__(self,name):
        self.fileName = name
    
        
class VtkUnstructured(OutputFile):
    
    def __init__(self,name):
        self.fileName = name
        self.ffmt = "%10.4f"
        
    def SetGrid(self,node,element):
        self.points = node
        self.cells = element
        
#    def AddNodeVector(self,d,dname):
        
#    def WriteScalarNodeData(self):
    
    def ElementType(self,e):
        
        if ( e.Type() == 'LINE2' ):
            return 3
        elif ( e.Type() == 'LINE3' ): 
            return 21  
        elif ( e.Type() == 'TRIA3' ): 
            return 5
        elif ( e.Type() == 'TRIA6' ): 
            return 22 
        elif ( e.Type() == 'QUAD4' ): 
            return 9
        elif ( e.Type() == 'QUAD8' ): 
            return 23
        elif ( e.Type() == 'TETRA4' ): 
            return 10
        elif ( e.Type() == 'TETRA10' ): 
            return 24
        elif ( e.Type() == 'HEXA8' ): 
            return 12
        elif ( e.Type() == 'HEXA20' ): 
            return 25
        
    def WriteFile(self):
        f = open(self.fileName,'w')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n\n')

        f.write('<UnstructuredGrid>\n')
        nn = len(self.points)
        ne = len(self.cells)
        f.write('<Piece NumberOfPoints="'+str(nn)+'" NumberOfCells="'+str(ne)+'">\n')
        
        # print point data 
        #   scalars
        #   vectors
        #   tensors
        
        # print cell data
        #   scalars
        #   vectors
        #   tensors
        
        # print points
        f.write('<Points>\n')
        f.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for n in self.points:
            for xi in n:
                f.write(self.ffmt % xi)
            for i in xrange(3-len(n)):
                f.write(self.ffmt % 0.0)
            f.write('\n')
        f.write('</DataArray>\n')
        f.write('</Points>\n')
        
        # print cells
        f.write('<Cells>\n')
        
        f.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for e in self.cells:
            for ni in e.conn:
                f.write(str(ni)+' ')
            f.write('\n')
        f.write('</DataArray>\n')
        
        
        f.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
        o=0
        c=0
        for e in self.cells:
            o=o+e.NumNodes()
            f.write(str(o)+' ')
            c=c+1
            if ( c%10 == 0 ):
                f.write('\n')
        if ( not c%10 == 0 ):
            f.write('\n')
        f.write('</DataArray>\n')
        
        f.write('<DataArray type="UInt8" Name="types" format="ascii">\n')        
        c=0
        for e in self.cells:
            o=o+e.NumNodes()
            f.write(str(self.ElementType(e))+' ')
            c=c+1
            if ( c%10 == 0 ):
                f.write('\n')
        if ( not c%10 == 0 ):
            f.write('\n')
        f.write('</DataArray>\n')
        
        f.write('</Cells>\n')
        
        f.write('</Piece>\n')
        f.write('</UnstructuredGrid>\n')
        f.write('</VTKFile>\n')
        f.close()
  
#-----------------------------------------------------------------------------              
import material as matl
import prop        
import meshing as mesh

plateWidth = 10.0
plateHeight = 10.0
he = .2

mat0 = matl.LinearElasticMat( 10e6, .33, 2.45e-4)
prop0 = prop.PlaneStress( mat0, thk = 1.0 )

nnx = int( plateWidth/he ) + 1
nny = int( plateHeight/he ) + 1
grid = mesh.MeshQuad4( np.array([[0,0],[plateWidth,0],[plateWidth,plateHeight],\
                 [0,plateHeight]], basic.FLOAT_TYPE ), nnx, nny, prop0 )
                 
vtkFile = VtkUnstructured('test3.vtu')
vtkFile.SetGrid( grid.node, grid.element )

vtkFile.WriteFile()