# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

# two elements same length, gap 2.0, always in contact to get the second solution

using FEMBase, MortarContact2DAD, Test, LinearAlgebra

X = Dict(
    1 => [0.0, 0.0],
    2 => [2.0, 0.0],
    3 => [0.0, 2.0],
    4 => [2.0, 2.0])

u = Dict(
    1 => [0.0, 0.0],
    2 => [0.0, 0.0],
    3 => [0.0, 0.0],
    4 => [0.0, 0.0])

function to_vec!(u, d)
    for (k, v) in d
        for i=1:2
            u[2*(k-1)+i] = v[i]
        end
    end
end

slave = Element(Seg2, (1, 2))
master = Element(Seg2, (3, 4))
update!((slave, master), "geometry", X)
update!((slave, master), "displacement", u)
problem = Problem(Contact2DAD, "test problem", 2, "displacement")
push!(problem.properties.always_in_contact, 1, 2)
add_slave_elements!(problem, slave)
add_master_elements!(problem, master)
problem.assembly.u = zeros(8)
problem.assembly.la = zeros(8)
to_vec!(problem.assembly.u, problem("displacement", 0.0))
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

eye = Matrix(1.0I, 4, 4)
K_expected = zeros(8, 8)
f_expected = zeros(8)
g_expected = zeros(8)
g_expected[1] = g_expected[3] = 2.0
C1_expected = zeros(8, 8)
C1_expected[1:4, :] .= [eye -eye]
C2_expected = zeros(8,8)
C2_expected[1:4, :] += [
 0.0  1.0  0.0  0.0  0.0 -1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0  0.0 -1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0]
C2_expected[1:4, :] += [
 1.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0]
D_expected = zeros(8,8)
D_expected[2,1] = D_expected[4,3] = 1.0

# @test isapprox(K, K_expected)
# @test isapprox(f, f_expected)
# @test isapprox(g, g_expected)
# @test isapprox(C1, C1_expected)
# @test isapprox(C2, C2_expected)
# @test isapprox(D, D_expected)

empty!(problem.assembly)
problem.properties.algorithm = 2
problem.properties.use_scaling = false
assemble!(problem, 0.0)

K, C1, C2, D, f, g = get_linear_system(problem, 8)
@test isapprox(K, K_expected)
@test isapprox(f, f_expected)
@test isapprox(g, g_expected)
@test isapprox(C1, C1_expected)
@test isapprox(D, D_expected)
@test isapprox(C2, C2_expected)
