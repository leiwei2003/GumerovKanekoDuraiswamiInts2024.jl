# two parallel triangles at various distances

using CompScienceMeshes
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test
using StaticArrays

# hypotenuses on a shared line

vts1 = [SVector(0.0, 1.0, 0.009), SVector(0.001, 1.0, 0.009), SVector(0.0, 1.0, 0.01)]
vts2 = [SVector(0.009, 0.0, 0.001), SVector(0.01, 0.0, 0.001), SVector(0.01, 0.0, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.989275221936278e-14]

vts1 = [SVector(0.0, 0.1, 0.009), SVector(0.001, 0.1, 0.009), SVector(0.0, 0.1, 0.01)]
vts2 = [SVector(0.009, 0.0, 0.001), SVector(0.01, 0.0, 0.001), SVector(0.01, 0.0, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.9734731539438897e-13]

vts1 = [SVector(0.0, 0.01, 0.009), SVector(0.001, 0.01, 0.009), SVector(0.0, 0.01, 0.01)]
vts2 = [SVector(0.009, 0.0, 0.001), SVector(0.01, 0.0, 0.001), SVector(0.01, 0.0, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.2287605621061589e-12]

vts1 = [SVector(0.0, 0.001, 0.009), SVector(0.001, 0.001, 0.009), SVector(0.0, 0.001, 0.01)]
vts2 = [SVector(0.009, 0.0, 0.001), SVector(0.01, 0.0, 0.001), SVector(0.01, 0.0, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.558500090581927e-12]

vts1 = [SVector(0.0, 0.0001, 0.009), SVector(0.001, 0.0001, 0.009), SVector(0.0, 0.0001, 0.01)]
vts2 = [SVector(0.009, 0.0, 0.001), SVector(0.01, 0.0, 0.001), SVector(0.01, 0.0, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.5632709331160492e-12]

# parallel hypotenuses

vts1 = [SVector(0.0, 0.0001, 0.0), SVector(0.001, 0.0001, 0.0), SVector(0.0, 0.0001, 0.001)]
vts2 = [SVector(0.009, 0.0, 0.01), SVector(0.01, 0.0, 0.01), SVector(0.01, 0.0, 0.009)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.5069241349342493e-12]

vts1 = [SVector(0.0, 0.001, 0.0), SVector(0.001, 0.001, 0.0), SVector(0.0, 0.001, 0.001)]
vts2 = [SVector(0.009, 0.0, 0.01), SVector(0.01, 0.0, 0.01), SVector(0.01, 0.0, 0.009)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.5026586342436529e-12]

vts1 = [SVector(0.0, 0.01, 0.0), SVector(0.001, 0.01, 0.0), SVector(0.0, 0.01, 0.001)]
vts2 = [SVector(0.009, 0.0, 0.01), SVector(0.01, 0.0, 0.01), SVector(0.01, 0.0, 0.009)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.2011159583783094e-12]

vts1 = [SVector(0.0, 0.1, 0.0), SVector(0.001, 0.1, 0.0), SVector(0.0, 0.1, 0.001)]
vts2 = [SVector(0.009, 0.0, 0.01), SVector(0.01, 0.0, 0.01), SVector(0.01, 0.0, 0.009)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.9723085272083167e-13]

vts1 = [SVector(0.0, 1.0, 0.0), SVector(0.001, 1.0, 0.0), SVector(0.0, 1.0, 0.001)]
vts2 = [SVector(0.009, 0.0, 0.01), SVector(0.01, 0.0, 0.01), SVector(0.01, 0.0, 0.009)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.9723085272083167e-13]