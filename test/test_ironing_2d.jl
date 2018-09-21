# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE.md

using FEMBase, MortarContact2DAD, Test, DelimitedFiles

using JuliaFEM

datadir = "test_ironing_2d"
meshfile = joinpath(datadir, "mesh1.med")
mesh = aster_read_mesh(meshfile)

create_node_set_from_element_set!(mesh, "TOOL_TO_BLOCK")
create_node_set_from_element_set!(mesh, "BLOCK_TO_TOOL")
S = collect(mesh.node_sets[:TOOL_TO_BLOCK])
M = collect(mesh.node_sets[:BLOCK_TO_TOOL])
N = setdiff(collect(keys(mesh.nodes)), [S;M])
renumber = Dict(j => i for (i,j) in enumerate([S;M;N]))
mesh.nodes = Dict(renumber[i] => j for (i,j) in mesh.nodes)
mesh.elements = Dict(i => [renumber[c] for c in j] for (i,j) in mesh.elements)
mesh.node_sets = Dict(n => Set(renumber[j] for j in v) for (n,v) in mesh.node_sets)
contact = Problem(Contact2DAD, "interface", 2, "displacement")
contact_slave_elements = create_elements(mesh, "TOOL_TO_BLOCK")
contact_master_elements = create_elements(mesh, "BLOCK_TO_TOOL")
add_slave_elements!(contact, contact_slave_elements)
add_master_elements!(contact, contact_master_elements)
contact.properties.rotate_normals = false
contact.properties.use_scaling = false
contact.properties.print_summary = true

get_mtx(name) = readdlm(joinpath(datadir, "ironing_2d_contact_$name.dat"))

contact.assembly.u = get_mtx("u")[:]
contact.assembly.la = get_mtx("la")[:]

function to_dict(u, ndofs, nnodes)
    return Dict(j => [u[ndofs*(j-1)+i] for i=1:ndofs] for j=1:nnodes)
end

nnodes = Int(length(contact.assembly.u)/2)
update!(contact_slave_elements, "displacement", 0.0 => to_dict(contact.assembly.u, 2, nnodes))
update!(contact_master_elements, "displacement", 0.0 => to_dict(contact.assembly.u, 2, nnodes))
update!(contact_slave_elements, "lambda", 0.0 => to_dict(contact.assembly.la, 2, nnodes))

contact.properties.iteration = 0
assemble!(contact, 0.0)

K = Matrix(contact.assembly.K)
K_expected = get_mtx("K")
@test isapprox(K, K_expected)
C1 = Matrix(contact.assembly.C1)
C1_expected = get_mtx("C1")
@test isapprox(C1, C1_expected)
C2 = Matrix(contact.assembly.C2)
C2_expected = get_mtx("C2")
@test isapprox(C2, C2_expected)
D = Matrix(contact.assembly.D)
D_expected = get_mtx("D")
@test isapprox(D, D_expected)
f = Vector(contact.assembly.f, 190)
f_expected = get_mtx("f")[:]
@test isapprox(f, f_expected; atol=1.0e-9)
g = Vector(contact.assembly.g, 30)
g_expected = get_mtx("g")[:]
@test isapprox(g, g_expected; atol=1.0e-9)
