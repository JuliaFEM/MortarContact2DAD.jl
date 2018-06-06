# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

using MortarContact2DAD
using Base.Test

@testset "MortarContact2DAD.jl tests" begin
    @testset "test_01.jl" begin include("test_01.jl") end
    @testset "test_02.jl" begin include("test_02.jl") end
    @testset "test_contact_1.jl" begin include("test_contact_1.jl") end
    @testset "test_contact_2.jl" begin include("test_contact_2.jl") end
end
