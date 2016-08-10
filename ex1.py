import numpy as np
import fempy

# This is a simple 2D plane stress example for illustration and debug purposes

# define the finite element mesh
pwidth=10
pheight=20
#mat0 = matl.LinearElasticMat( 10e6, .33, 2.45e-4 )
#prop0 = prop.PlaneStress( mat0, 1.0 )
#fegrid = mesh.MeshQuad4( np.array([[0,0],[pwidth,0],[pwidth,pheight],[0,pheight]]), 3, 5, prop0 )
#
## setup a simple dofmap with 2 dof per node
#dofMap = dmap.FixedDofMap(2)
#formulation = prob.ProblemPlaneStress()
#
## compute the stiffness matrix
#Kmat = formulation.ComputeStiffness( fegrid.node, fegrid.element, dofMap )
#
## define the boundary conditions
#spcs = bcs.EssentialBCs()
#spcs.AddSpcs( fegrid.EdgeNIDs(0), [1], dofmap=dofMap )  # fix bottom edge in y
#spcs.AddSpcs( fegrid.EdgeNIDs(3), [0], dofmap=dofMap )  # fix left edge in x
#
#loads = bcs.NaturalBCs()
#loads.AddTraction( fegrid.EdgeElem(2), [0.0, 100.0] )
#
## solve the system
#solver = basic.FeSolver(Kmat)
#d, freac = solver.Solve( loads.Fext(fegrid.node,dofMap), spcs.SPCs() )
#
## compute stress and write results
#stress,svm,strain = formulation.ComputeStress( fegrid.node, fegrid.element, d, dofMap )
#fedata = out.FeaData(fegrid.node,fegrid.element)
#fedata.SetData(d,dofMap,[0,1],stress,svm)
#fedata.PlotScalar()
#fedata.WriteVtkFile('plot.vtu')