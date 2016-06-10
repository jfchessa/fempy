#!/usr/bin/python
#
# A 3D continuum structrual finite element code

import numpy as np
import basic
import dofmap as dmap
import material as matl
import prop
import meshing as mesh
import linear_problems as prob
import bcs
import output as out

#import profile
#profile.run('execfile("ex2.py"); print')


print '********************* 3D FINITE ELEMENT PROGRAM *********************'

print 'Initializing grid'
# define the finite element mesh
pwidth=20
plength=20
pheight=5
mat0 = matl.LinearElasticMat( 10e6, .33, 2.45e-4 )  # linear elastic isotropic material
prop0 = prop.Solid3D( mat0 )                        # 3D solid property
fegrid = mesh.MeshHexa8( np.array([ [0,0,0], [plength,0,0], [plength,pwidth,0], [0,pwidth,0], \
        [0,0,pheight], [plength,0,pheight], [plength,pwidth,pheight], [0,pwidth,pheight] ]),\
        21, 21, 6, prop0 )

# setup a simple dofmap with 3 dof per node
dofMap = dmap.FixedDofMap(3)
formulation = prob.Problem3D()

print 'Computing stiffness matrix'
# compute the stiffness matrix
Kmat = formulation.ComputeStiffness( fegrid.node, fegrid.element, dofMap )

# define symmetry boundary conditions
spcs = bcs.EssentialBCs()
spcs.AddSpcs( fegrid.FaceNIDs(1), [0], dofmap=dofMap )  # fix -x face in x
spcs.AddSpcs( fegrid.FaceNIDs(3), [1], dofmap=dofMap )  # fix -y face in y
spcs.AddSpcs( fegrid.FaceNIDs(5), [2], dofmap=dofMap )  # fix -z face in z

# define rhs load
loads = bcs.NaturalBCs()
loads.AddTraction( fegrid.FaceElem(0), [100.0, 0.0, 0.0] )  # load on the +x face 
rhs = loads.Fext(fegrid.node,dofMap)

#print 'Solving Kd=f'
# solve the system
solver = basic.FeSolver(Kmat)
d, freac = solver.Solve( rhs, spcs.SPCs() )

print 'Writting results'
# compute stress and write results
stress,svm,strain = formulation.ComputeStress( fegrid.node, fegrid.element, d, dofMap )
fedata = out.FeaData(fegrid.node,fegrid.element)
fedata.SetData(d,dofMap,[0,1,2],stress,svm)
fedata.WriteVtkFile('plot.vtu')

print '**************************** END OF RUN *****************************'