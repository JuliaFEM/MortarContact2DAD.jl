# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

"""
    2D mortar contact mechanics for JuliaFEM using Automatic Differentiation
"""
module MortarContact2DAD

using ForwardDiff

using Reexport
@reexport using FEMBase
import FEMBase: get_unknown_field_name, assemble_elements!, add_elements!

const MortarElements2D = Union{Seg2,Seg3}

include("mortar2dad.jl")
include("contact2dad.jl")
export Mortar2DAD, Contact2DAD
export add_slave_elements!, add_master_elements!
export get_slave_elements, get_master_elements
export get_slave_dofs, get_master_dofs

end
