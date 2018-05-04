# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

"""
    2D mortar contact mechanics for JuliaFEM using Automatic Differentiation
"""
module MortarContact2DAD

using Reexport
@reexport using FEMBase

export Contact2DAD
export add_slave_elements!, add_master_elements!
export get_slave_elements, get_master_elements
export get_slave_dofs, get_master_dofs

end
