# test that are currecntly still problematic

using CompScienceMeshes
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test
using StaticArrays

vts1 = [SVector(0.001, 0.01, 0.0), SVector(0.002, 0.01, 0.0), SVector(0.001, 0.01, 0.001)]
vts2 = [SVector(0.0, 0.001, 0.0), SVector(0.001, 0.0, 0.0), SVector(0.001, 0.001, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh1,mesh2) ≈ [2.1239384709441043e-12]

vts1 = [SVector(0.0, 0.0, 0.006), SVector(0.001, 0.0, 0.006), SVector(0.001, 0.0, 0.005)]
vts2 = [ SVector(0.0, 0.001, 0.0),SVector(0.001, 0.0, 0.0),SVector(0.001, 0.001, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [3.4839699355741572e-12]

vts1 = [SVector(0.002, 0.01, 0.001), SVector(0.003, 0.01, 0.0), SVector(0.003, 0.01, 0.001)]
vts2 = [SVector(0.01, 0.01, 0.0),SVector(0.01, 0.009000000000000001, 0.001),SVector(0.01, 0.01, 0.001)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [2.7085862026852206e-12]

# parallel triangles

vts1 = 10* [SVector(0.0, 0.01, 0.009000000000000001), SVector(0.001, 0.01, 0.009000000000000001), SVector(0.0, 0.01, 0.01)]
vts2 = [SVector(0.009000000000000001, 0.0, 0.001), SVector(0.01, 0.0, 0.001), SVector(0.01, 0.0, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh2,mesh1) ≈ [1.2287605621063186e-12]