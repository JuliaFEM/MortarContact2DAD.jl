# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE.md

using JuliaFEM, MortarContact2DAD

#datadir = Pkg.dir("MortarContact2DAD", "examples", "2d_ironing")
datadir = abspath(joinpath(dirname(pathof(MortarContact2DAD)), "..", "examples", "2d_ironing"))
meshfile = joinpath(datadir, "mesh0.med")

mesh = aster_read_mesh(meshfile)

# nel = length(mesh.element_sets[:BLOCK_TO_TOOL])
# @info("Number of elements in set BLOCK_TO_TOOL: $nel")
# for elid in mesh.element_sets[:BLOCK_TO_TOOL]
#     for nid in mesh.elements[elid]
#         x, y = mesh.nodes[nid]
#         if x > 2.2
#             setdiff!(mesh.element_sets[:BLOCK_TO_TOOL], elid)
#         end
#     end
# end
# nel = length(mesh.element_sets[:BLOCK_TO_TOOL])
# @info("Number of elements in set BLOCK_TO_TOOL after filter: $nel")

create_node_set_from_element_set!(mesh, "TOOL_TO_BLOCK")
create_node_set_from_element_set!(mesh, "BLOCK_TO_TOOL")
S = collect(mesh.node_sets[:TOOL_TO_BLOCK])
M = collect(mesh.node_sets[:BLOCK_TO_TOOL])
N = setdiff(collect(keys(mesh.nodes)), [S;M])
renumber = Dict(j => i for (i,j) in enumerate([S;M;N]))
mesh.nodes = Dict(renumber[i] => j for (i,j) in mesh.nodes)
mesh.elements = Dict(i => [renumber[c] for c in j] for (i,j) in mesh.elements)
mesh.node_sets = Dict(n => Set(renumber[j] for j in v) for (n,v) in mesh.node_sets)
S = collect(mesh.node_sets[:TOOL_TO_BLOCK])
M = collect(mesh.node_sets[:BLOCK_TO_TOOL])
N = setdiff(collect(keys(mesh.nodes)), [S;M])

tool = Problem(mesh, Elasticity, "TOOL", 2)
tool.properties.formulation = :plane_stress
tool.properties.finite_strain = true
tool.properties.geometric_stiffness = true
update!(tool, "youngs modulus", 6896.0)
update!(tool, "poissons ratio", 0.32)

block = Problem(mesh, Elasticity, "BLOCK", 2)
block.properties.formulation = :plane_stress
block.properties.finite_strain = true
block.properties.geometric_stiffness = true
update!(block, "youngs modulus", 689.6)
update!(block, "poissons ratio", 0.32)

bc1 = Problem(mesh, Dirichlet, "FIXED", 2, "displacement")
update!(bc1, "displacement 1", 0.0)
update!(bc1, "displacement 2", 0.0)

bc2 = Problem(mesh, Dirichlet, "TOP", 2, "displacement")
update!(bc2, "displacement 1", 0.0 => 0.0)
update!(bc2, "displacement 2", 0.0 => 0.0)
update!(bc2, "displacement 1", 1.8 => 0.0)
update!(bc2, "displacement 2", 1.8 => -1.0)
update!(bc2, "displacement 1", 11.4 => 9.60)
update!(bc2, "displacement 2", 11.4 => -1.0)
update!(bc2, "displacement 1", 13.2 => 9.60)
update!(bc2, "displacement 2", 13.2 => 0.00)

contact = Problem(Contact2DAD, "interface", 2, "displacement")
contact_slave_elements = create_elements(mesh, "TOOL_TO_BLOCK")
contact_master_elements = create_elements(mesh, "BLOCK_TO_TOOL")
add_slave_elements!(contact, contact_slave_elements)
add_master_elements!(contact, contact_master_elements)
contact.properties.rotate_normals = true
contact.properties.use_scaling = true
contact.properties.print_summary = true
contact.assembly.u = zeros(2*length(mesh.nodes))
contact.assembly.la = zeros(2*length(mesh.nodes))

analysis = Analysis(Nonlinear)
analysis.properties.max_iterations = 30
add_problems!(analysis, [tool, block, bc1, bc2, contact])
xdmf = Xdmf(joinpath(datadir, "2d_ironing_100"); overwrite=true)
add_results_writer!(analysis, xdmf)

end_time = 13.2
timespan = range(0.0, stop=end_time, length=Int(end_time*60))
for analysis.properties.time in timespan
    contact.properties.iteration = 0
    run!(analysis)
end
close(xdmf.hdf)

using Plots

xm = [1.86784, 3.99604]
xs1 = [0.869307, 3.65322]
xs2 = [1.07536, 3.55731]
ns1 = [-0.413177, -0.910651]
ns2 = [-0.172203, -0.985061]
plt = scatter!([xm[1]], [xm[2]], markersize=6, c=:red, legend=false)
plot!(plt, [xs1[1], xs2[1]], [xs1[2], xs2[2]], c=:blue, w=2, m=:dot)
plot!(plt, [xs1[1], xs1[1]+ns1[1]], [xs1[2], xs1[2]+ns1[2]], c=:black, w=1)
plot!(plt, [xs2[1], xs2[1]+ns2[1]], [xs2[2], xs2[2]+ns2[2]], c=:black, w=1)

tim = first(contact_slave_elements)("displacement").data[end-1].first
#tim = analysis.properties.time

function plot_surfaces(time)

    for element in contact_slave_elements
        (X1, Y1), (X2, Y2) = element("geometry", time)
        (u1, v1), (u2, v2) = element("displacement", time)
        x1 = X1 + u1
        y1 = Y1 + v1
        x2 = X2 + u2
        y2 = Y2 + v2
        plot!(plt, [x1, x2], [y1, y2], c=:blue, w=2, m=:dot)
    end

    for element in contact_master_elements
        (X1, Y1), (X2, Y2) = element("geometry", time)
        (u1, v1), (u2, v2) = element("displacement", time)
        #u1, v1, u2, v2 = contact.assembly.u[get_gdofs(contact, element)]
        x1 = X1 + u1
        y1 = Y1 + v1
        x2 = X2 + u2
        y2 = Y2 + v2
        plot!(plt, [x1, x2], [y1, y2], c=:red, w=2, m=:dot)
    end

end

plot_surfaces(tim)
display(plt)
