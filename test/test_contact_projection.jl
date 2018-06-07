# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

using MortarContact2DAD
using MortarContact2DAD: project_from_master_to_slave_ad
using ForwardDiff
using FEMBase.Test

function func1(u)
    x1 = [0.0, 0.0] + u[1:2]
    x2 = [0.0, 2.0] + u[3:4]
    x3 = [0.5, 0.5] + u[5:6]
    t1_ = (x2-x1) / norm(x2-x1)
    n1 = n2 = [0.0 1.0; -1.0 0.0] * t1_
    x1_ = DVTI(x1, x2)
    n1_ = DVTI(n1, n2)
    slave = Element(Seg2, [1, 2])
    xi1 = project_from_master_to_slave_ad(slave, x1_, n1_, x3; debug=true)
    return [xi1]
end

J = ForwardDiff.jacobian(func1, zeros(6))
@test isapprox(J, [-0.25 -0.75 0.25 -0.25 0.0 1.0])

function func2(u)
    x1 = [0.0, 0.0] + u[1:2]
    x2 = [0.0, 2.0] + u[3:4]
    x3 = [0.5, 0.5] + u[5:6]
    t1_ = (x2-x1) / norm(x2-x1)
    n1 = n2 = [0.0 1.0; -1.0 0.0] * t1_
    x1_ = DVTI(x1, x2)
    n1_ = DVTI(n1, n2)
    slave = Element(Seg2, [1, 2])
    xi1 = project_from_master_to_slave_ad(slave, x1_, n1_, x3; debug=true, max_iterations=0)
    return [xi1]
end

@test_throws Exception ForwardDiff.jacobian(func2, zeros(6))
