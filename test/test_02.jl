# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

using MortarContact2DAD
using FEMBase.Test
import FEMBase.Test: @test_resource

X = Dict(1 => [-2.0, 0.0], 2 => [ 0.0, 0.0], 3 => [ 2.0, 0.0],
         4 => [-2.0, 0.0], 5 => [ 0.0, 0.0], 6 => [-2.0, -2.0],
         7 => [0.0, -2.0], 8 => [2.0, -2.0], 9 => [-2.0, 2.0],
         10 => [0.0, 2.0])

 D = [0.5, 0.0]

for j in [4, 5, 9, 10]
    X[j] += D
end

contact = Problem(Contact2DAD, "contact", 2, "displacement")
element7 = Element(Seg2, [1, 2])
element8 = Element(Seg2, [2, 3])
element9 = Element(Seg2, [4, 5])
contact_slave_elements = [element7, element8]
contact_master_elements = [element9]
add_slave_elements!(contact, contact_slave_elements)
add_master_elements!(contact, contact_master_elements)
all_elements = [element7, element8, element9]
update!(all_elements, "geometry", X)

using JLD
i1 = load(@test_resource("iter_1.jld"))
i2 = load(@test_resource("iter_2.jld"))
i3 = load(@test_resource("iter_3.jld"))

function to_dict(u, ndofs, nnodes)
    return Dict(j => [u[ndofs*(j-1)+i] for i=1:ndofs] for j=1:nnodes)
end

ndofs = 20
contact.properties.iteration = 1
for (data, time) in zip([i1, i2, i3], [0.0, 1.0, 2.0])
    println("Testing data for iteration / time $time ")
    empty!(contact.assembly)
    contact.assembly.u = data["u"]
    contact.assembly.la = data["la"]
    update!(contact, "displacement", time => to_dict(data["u"], 2, 5))
    update!(contact, "lambda", time => to_dict(data["la"], 2, 5))
    assemble!(contact, time)
    K = full(contact.assembly.K, ndofs, ndofs)
    C1 = full(contact.assembly.C1, ndofs, ndofs)
    C2 = full(contact.assembly.C2, ndofs, ndofs)
    D = full(contact.assembly.D, ndofs, ndofs)
    f = full(contact.assembly.f, ndofs, 1)
    g = full(contact.assembly.g, ndofs, 1)
    @test isapprox(K, data["K"])
    @test isapprox(C1, data["C1"])
    @test isapprox(C2, data["C2"])
    @test isapprox(D, data["D"])
    @test isapprox(f, data["f"])
    @test isapprox(g, data["g"])
end
