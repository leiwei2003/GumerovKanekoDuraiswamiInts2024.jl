# tests for the double layer operator that are currently still problematic

using CompScienceMeshes
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test
using StaticArrays

vts1 = [SVector(0.001, 0.01, 0.0), SVector(0.002, 0.01, 0.0), SVector(0.001, 0.01, 0.001)]
vts2 = [SVector(0.0, 0.001, 0.0), SVector(0.001, 0.0, 0.0), SVector(0.001, 0.001, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh1,mesh2,operator="doublelayer") ≈ [ 8.081884922746134e-12]

vts1 = [SVector(0.0, 0.0, 0.006), SVector(0.001, 0.0, 0.006), SVector(0.001, 0.0, 0.005)]
vts2 = [ SVector(0.0, 0.001, 0.0),SVector(0.001, 0.0, 0.0),SVector(0.001, 0.001, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh1,mesh2,operator="doublelayer") ≈ [6.055001550474312e-10]

vts1 = [SVector(0.002, 0.01, 0.001), SVector(0.003, 0.01, 0.0), SVector(0.003, 0.01, 0.001)]
vts2 = [SVector(0.01, 0.01, 0.0),SVector(0.01, 0.009000000000000001, 0.001),SVector(0.01, 0.01, 0.001)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh1,mesh2,operator="doublelayer") ≈ [3.681758955775994e-10]

vts1 = [SVector(0.0, 0.01, 0.009000000000000001), SVector(0.001, 0.01, 0.009000000000000001), SVector(0.0, 0.01, 0.01)]
vts2 = [SVector(0.009000000000000001, 0.0, 0.001), SVector(0.01, 0.0, 0.001), SVector(0.01, 0.0, 0.0)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh1,mesh2,operator="doublelayer") ≈ [4.693035723757102e-11]

vts1 = [SVector(1.0, -1.0, -0.6666666666657435), SVector(1.0, -0.7223240424047911, -0.5045874500029877), SVector(1.0, -1.0, -0.3333333333314863)]
vts2 = [SVector(-1.0, -0.7261013160371913, -0.5057662089637784), SVector(-1.0, -0.7559830641446982, -0.755983064144203), SVector(-1.0, -1.0, -0.6666666666657435)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh1,mesh2,operator="doublelayer") ≈ [ -2.8885991615515054e-5]

vts1 = [SVector(-0.7219173941427399, 1.0, -0.4911910680760563), SVector(-0.7559830641446982, 1.0, -0.755983064144203), SVector(-1.0, 1.0, -0.6666666666657435)]
vts2 = [SVector(-0.7223240424047911, -1.0, -0.5045874500029878), SVector(-1.0, -1.0, -0.6666666666657435), SVector(-0.7559830641430324, -1.0, -0.7559830641430324)]
mesh1 = Mesh(vts1, [SVector(1,2,3)])
mesh2 = Mesh(vts2, [SVector(1,2,3)])
@test IntegrateMesh(mesh1,mesh2,operator="doublelayer") ≈ [-2.154576551327468e-5]