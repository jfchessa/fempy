import numpy as np
import scipy as sp
from scipy import sparse

import basic
import prop
import material as matl
import element as elem 
import bcs
import meshing as mesh  
import output as plot
import linear_problems as prob

import pdb
        
length =10.0
width = 10.0
height = 2.0
he=0.2
mat0 = matl.LinearElasticMat( 10e6, .33, 2.45e-4)
prop0 = prop.Solid3D( mat0 )

nnx = int( length/he ) + 1
nny = int( width/he ) + 1
nnz = int( height/he ) + 1

corners = np.array([[0,0,0],[length,0,0],[length,width,0],[0,width,0],\
                [0,0,height],[length,0,height],[length,width,height],\
                [0,width,height]], basic.FLOAT_TYPE )

grid = mesh.MeshHexa8( corners, nnx, nny, nnz, prop0 )
 
# compute the stiffness matrix                
problem = prob.Problem3D()
Kmat = problem.ComputeStiffness( grid.node, grid.element )

# compute the external force
rightFace = grid.FaceElem(0)
trac = np.array([ 1000.0, 0.0, 0.0])
fext = bcs.compute_fext_tract(grid.node,rightFace,trac)

# set the spcs
ifix = np.concatenate( (3*grid.FaceNIDs(1), 3*grid.FaceNIDs(3)+1, 3*grid.FaceNIDs(5)+2) )

# solve Kd=f
d, freac = basic.fesolve(Kmat,fext,ifix) 

# compute stress
stress,svm,strain = problem.ComputeStress( grid.node, grid.element, d )

# plot the results
plotData = plot.FeaData(grid.node,grid.element)
plotData.SetDisplacement(d)
plotData.SetStress(stress)
plotData.SetEffectiveStress(svm)
plotData.PlotScalar()
plotData.WriteVtkFile('plot.vtu')