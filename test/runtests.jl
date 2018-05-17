# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

using MortarContact2DAD
using Base.Test

@testset "MortarContact2DAD.jl tests" begin
    @testset "test_01.jl" begin include("test_01.jl") end
end
