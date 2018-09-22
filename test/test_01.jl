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
problem = Problem(Contact2DAD, "test problem", 2, "displacement")
add_slave_elements!(problem, slave)
add_master_elements!(problem, master)
problem.assembly.u = zeros(8)
problem.assembly.la = zeros(8)
assemble!(problem, 0.0)

function get_linear_system(problem, ndofs)
    K = Matrix(problem.assembly.K, ndofs, ndofs)
    C1 = Matrix(problem.assembly.C1, ndofs, ndofs)
    C2 = Matrix(problem.assembly.C2, ndofs, ndofs)
    D = Matrix(problem.assembly.D, ndofs, ndofs)
    f = Vector(problem.assembly.f, ndofs)
    g = Vector(problem.assembly.g, ndofs)
    return K, C1, C2, D, f, g
end

K, C1, C2, D, f, g = get_linear_system(problem, 8)

@test isapprox(K, zeros(8, 8))
@test isapprox(f, zeros(8))
@test isapprox(g, zeros(8))

C1_expected = [
 1.0  0.0  0.0  0.0  -1.0   0.0   0.0   0.0
 0.0  1.0  0.0  0.0   0.0  -1.0   0.0   0.0
 0.0  0.0  1.0  0.0   0.0   0.0  -1.0   0.0
 0.0  0.0  0.0  1.0   0.0   0.0   0.0  -1.0]
@test isapprox(C1[1:4,:], C1_expected)

C2_expected = [
 0.0  1.0  0.0  0.0  0.0  -1.0  0.0   0.0
 0.0  0.0  0.0  0.0  0.0   0.0  0.0   0.0
 0.0  0.0  0.0  1.0  0.0   0.0  0.0  -1.0
 0.0  0.0  0.0  0.0  0.0   0.0  0.0   0.0]
@test isapprox(C2[1:4,:], C2_expected)

D_expected = [
 0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0]
@test isapprox(D[1:4,1:4], D_expected)

empty!(problem.assembly)
problem.properties.algorithm = 2
problem.properties.use_scaling = false
assemble!(problem, 0.0)

K, C1, C2, D, f, g = get_linear_system(problem, 8)
@test isapprox(K, zeros(8, 8))
@test isapprox(f, zeros(8))
@test isapprox(g, zeros(8))
@test isapprox(C1[1:4,:], C1_expected)
@test isapprox(C2[1:4,:], C2_expected)
@test isapprox(D[1:4,1:4], D_expected)
