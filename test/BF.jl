# two parallel triangles at various distances

using CompScienceMeshes
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test
using StaticArrays

#BigFloat shortcut
BF = BigFloat

# hypotenuses on a shared line

vts1 = [SVector(BF(0.0), BF(100), BF(0.009)), SVector(BF(0.001), BF(100), BF(0.009)), SVector(BF(0.0), BF(100), BF(0.01))]
vts2 = [SVector(BF(0.009), BF(0.0), BF(0.001)), SVector(BF(0.01), BF(0.0), BF(0.001)), SVector(BF(0.01), BF(0.0), BF(0.0))]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.9894367724900594e-16]

vts1 = [SVector(0.0, 100, 0.009), SVector(0.001, 100, 0.009), SVector(0.0, 100, 0.01)]
vts2 = [SVector(0.009, 0.0, 0.001), SVector(0.01, 0.0, 0.001), SVector(0.01, 0.0, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.9894367724900594e-16]
