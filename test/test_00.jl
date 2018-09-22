# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

# Fundamentals

using FEMBase, MortarContact2DAD, Test
using MortarContact2DAD: get_slave_dofs, get_master_dofs
using MortarContact2DAD: project_from_master_to_slave_ad, project_from_slave_to_master_ad

X = Dict(
    1 => [0.0, 0.0],
    2 => [1.0, 0.0],
    3 => [0.0, 1.0],
    4 => [1.0, 1.0])

u = Dict(
    1 => [0.0, 0.0],
    2 => [0.0, 0.0],
    3 => [0.0, 0.0],
    4 => [0.0, 0.0])

slave = Element(Seg2, (1, 2))
master = Element(Seg2, (3, 4))
update!((slave, master), "geometry", X)
update!((slave, master), "displacement", u)
update!(slave, "normal", ([0.0, 1.0], [0.0, 1.0]))
problem = Problem(Mortar2DAD, "test problem", 2, "displacement")

@test_throws Exception add_elements!(problem, [slave])
add_slave_elements!(problem, [slave])
add_master_elements!(problem, [master])
@test get_slave_dofs(problem) == [1, 2, 3, 4]
@test get_master_dofs(problem) == [5, 6, 7, 8]
