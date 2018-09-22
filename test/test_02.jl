# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

# two elements same length, gap 0

using FEMBase, MortarContact2DAD, Test

X = Dict(
    1 => [0.0, 0.0],
    2 => [2.0, 0.0],
    3 => [0.0, 0.0],
    4 => [2.0, 0.0])

u = Dict(
    1 => [0.0, 0.0],
    2 => [0.0, 0.0],
    3 => [0.0, 0.0],
    4 => [0.0, 0.0])

slave = Element(Seg2, (1, 2))
master = Element(Seg2, (3, 4))
update!((slave, master), "geometry", X)
update!((slave, master), "displacement", u)
problem = Problem(Mortar2DAD, "test problem", 2, "displacement")
add_slave_elements!(problem, slave)
add_master_elements!(problem, master)
problem.assembly.u = zeros(8)
problem.assembly.la = zeros(8)
assemble!(problem, 0.0)
C1_expected = 1/6*[
               2.0 0.0 1.0 0.0 -2.0 0.0 -1.0 0.0
               0.0 2.0 0.0 1.0 0.0 -2.0 0.0 -1.0
               1.0 0.0 2.0 0.0 -1.0 0.0 -2.0 0.0
               0.0 1.0 0.0 2.0 0.0 -1.0 0.0 -2.0]
C1 = problem.assembly.C1
C2 = problem.assembly.C2
@test isapprox(C1, C2)
@test isapprox(C1, C1_expected)
