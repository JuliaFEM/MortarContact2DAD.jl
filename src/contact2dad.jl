# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

using ForwardDiff: value

mutable struct Contact2DAD <: BoundaryProblem
    master_elements :: Vector{Element}
    rotate_normals :: Bool
    dual_basis :: Bool
    iteration :: Int
    contact_state_in_first_iteration :: Symbol
    always_in_contact :: Set{Int64}
    use_scaling :: Bool
    alpha :: Dict{Int64, Real}
    beta :: Dict{Int64, Real}
    nadj_nodes :: Dict{Int64, Real}
    scaling_factors :: Dict{Int64, Real}
    print_summary :: Bool
    contact_pairs :: Dict{Element, Set{Element}}
    algorithm :: Int
end

function Contact2DAD()
    master_elements = []
    rotate_normals = false
    dual_basis = true
    iteration = 0
    contact_state_in_first_iteration = :AUTO
    always_in_contact = Set()
    use_scaling = true
    alpha = Dict()
    beta = Dict()
    nadj_nodes = Dict()
    scaling_factors = Dict()
    print_summary = true
    contact_pairs = Dict()
    algorithm = 1
    return Contact2DAD(master_elements, rotate_normals, dual_basis, iteration,
                       contact_state_in_first_iteration, always_in_contact,
                       use_scaling, alpha, beta, nadj_nodes, scaling_factors,
                       print_summary, contact_pairs, algorithm)
end

function FEMBase.add_elements!(::Problem{Contact2DAD}, ::Any)
    error("use `add_slave_elements!` and `add_master_elements!` to add ",
          "elements to the Mortar2D problem.")
end

function FEMBase.add_slave_elements!(problem::Problem{Contact2DAD}, elements...)
    for element in elements
        push!(problem.elements, element)
    end
end

function FEMBase.add_master_elements!(problem::Problem{Contact2DAD}, elements...)
    for element in elements
        push!(problem.properties.master_elements, element)
    end
end

function FEMBase.add_slave_elements!(problem::Problem{Contact2DAD}, element_sets::Union{Vector, Tuple}...)
    for element_set in element_sets
        add_slave_elements!(problem, element_set...)
    end
end

function FEMBase.add_master_elements!(problem::Problem{Contact2DAD}, element_sets::Union{Vector, Tuple}...)
    for element_set in element_sets
        add_master_elements!(problem, element_set...)
    end
end

function FEMBase.get_slave_elements(problem::Problem{Contact2DAD})
    return problem.elements
end

function FEMBase.get_master_elements(problem::Problem{Contact2DAD})
    return problem.properties.master_elements
end

function FEMBase.get_elements(problem::Problem{Contact2DAD})
    return [get_slave_elements(problem); get_master_elements(problem)]
end

function FEMBase.get_master_elements(problem::Problem{Contact2DAD}, slave_element::Element)
    return problem.properties.contact_pairs[slave_element]
end

function get_dofs(problem::Problem{Contact2DAD}, elements)
    dofs = Int64[]
    for element in elements
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

function get_slave_dofs(problem::Problem{Contact2DAD})
    return get_dofs(problem, get_slave_elements(problem))
end

function get_master_dofs(problem::Problem{Contact2DAD})
    return get_dofs(problem, get_master_elements(problem))
end

function get_nodes(elements)
    nodes = Set{Int64}()
    for element in elements
        for j in get_connectivity(element)
            push!(nodes, j)
        end
    end
    return nodes
end

function get_slave_nodes(problem::Problem{Contact2DAD})
    return get_nodes(get_slave_elements(problem))
end

function get_master_nodes(problem::Problem{Contact2DAD})
    return get_nodes(get_master_elements(problem))
end

function project_from_slave_to_master(::Type{Val{:Seg2}}, x1, n1, xm1, xm2)
    x11, x12 = x1
    n11, n12 = n1
    x31, x32 = xm1
    x41, x42 = xm2
    nom = -2*n11*x12 + n11*x32 + n11*x42 + 2*n12*x11 - n12*x31 - n12*x41
    denom = n11*x32 - n11*x42 - n12*x31 + n12*x41
    return nom/denom
end

function project_from_master_to_slave(::Type{Val{:Seg2}}, xm, xs1, xs2, ns1, ns2)

    x11, x12 = xs1
    x21, x22 = xs2
    x31, x32 = xm
    n11, n12 = ns1
    n21, n22 = ns2
    a = (-n11*x12 + n11*x22 + n12*x11 - n12*x21 + n21*x12 - n21*x22 - n22*x11 + n22*x21)/4
    if isapprox(a, 0.0; atol=1.0e-3)
        nom = n11*x12 + n11*x22 - 2*n11*x32 - n12*x11 - n12*x21 + 2*n12*x31
        denom = n11*x12 - n11*x22 - n12*x11 + n12*x21
        return nom/denom
    end
    b = -(-n11*x12 + n11*x32 + n12*x11 - n12*x31 + n21*x22 - n21*x32 - n22*x21 + n22*x31)/2
    c = (-n11*x12 - n11*x22 + 2*n11*x32 + n12*x11 + n12*x21 - 2*n12*x31 - n21*x12 - n21*x22 + 2*n21*x32 + n22*x11 + n22*x21 - 2*n22*x31)/4
    d = b^2 - 4*a*c
    if d < 0.0
        xm = map(value, xm)
        xs1 = map(value, xs1)
        xs2 = map(value, xs2)
        ns1 = map(value, ns1)
        ns2 = map(value, ns2)
        a = value(a)
        b = value(b)
        c = value(c)
        d = value(d)
        @error("xm = $xm")
        @error("xs1 = $xs1")
        @error("xs2 = $xs2")
        @error("ns1 = $ns1")
        @error("ns2 = $ns2")
        @error("a = $a, b = $b, c = $c, d = $d")
        error("When projecting from master to slave, got error. Negative discriminant.")
    end
    sol1 = 1.0/(2.0*a)*(-b+sqrt(d))
    sol2 = 1.0/(2.0*a)*(-b-sqrt(d))
    if abs(sol1) < abs(sol2)
        return sol1
    else
        return sol2
    end
end

const quadrature_points = ((0.23692688505618908, -0.906179845938664),
                           (0.47862867049936647, -0.538469310105683),
                           (0.56888888888888890,  0.000000000000000),
                           (0.47862867049936647,  0.538469310105683),
                           (0.23692688505618908,  0.906179845938664))

function calculate_dual_basis_coefficients!(Ae::Matrix{T}, problem, slave_element, x, n) where {T}

    cs1, cs2 = get_connectivity(slave_element)

    De = zeros(T, 2, 2)
    Me = zeros(T, 2, 2)
    nsegments = 0

    for master_element in get_master_elements(problem, slave_element)

        cm1, cm2 = get_connectivity(master_element)

        # calculate segmentation: we care only about endpoints
        xi1a = project_from_master_to_slave(Val{:Seg2}, x[cm1], x[cs1], x[cs2], n[cs1], n[cs2])
        xi1b = project_from_master_to_slave(Val{:Seg2}, x[cm2], x[cs1], x[cs2], n[cs1], n[cs2])
        xi1a = clamp(xi1a, -1.0, 1.0)
        xi1b = clamp(xi1b, -1.0, 1.0)
        l = 1/2*abs(xi1b-xi1a)
        isapprox(l, 0.0) && continue # no contribution in this master element

        nsegments += 1
        for (weight, xi) in quadrature_points
            detj = norm(0.5*(x[cs2]-x[cs1]))
            w = weight*detj*l
            xi_s = 0.5*(1.0-xi)*xi1a + 0.5*(1.0+xi)*xi1b
            N1 = [0.5*(1.0-xi_s), 0.5*(1.0+xi_s)]
            De += w*Matrix(Diagonal(N1))
            Me += w*kron(N1, N1')
        end
    end

    if nsegments == 0
        Ae[:,:] .= Matrix(1.0I, 3, 3)
        return nothing
    else
        a, b, c, d = Me
        invMe = 1.0/(a*d-b*c) * [d -c; -b a]
        Ae[:,:] .= De*invMe
        return nothing
    end

end

function assemble_elements_preprocess!(problem::Problem{Contact2DAD}, assembly::Assembly,
                                       elements::Vector{Element{Seg2}}, time::Float64)

    @debug("First iteration, creating contact pairing in preprocess.")

    S = get_slave_nodes(problem)
    slave_elements = get_slave_elements(problem)
    props = problem.properties

    @debug("Calculate nodal normals")
    normals = Dict{Int64, Vector{Float64}}(nid => zeros(Float64, 2) for nid in S)
    coeff = props.rotate_normals ? -1.0 : 1.0
    for slave_element in slave_elements
        conn = get_connectivity(slave_element)
        X1 = slave_element("geometry", time)
        u1 = slave_element("displacement", time)
        # @info("slave element info", X1, u1)
        x1, x2 = tuple((Xi+ui for (Xi,ui) in zip(X1,u1))...)
        t1, t2 = coeff * (x2-x1) / norm(x2-x1)
        n = [-t2, t1]
        for c in conn
            normals[c] += n
        end
    end

    for nid in keys(normals)
        normals[nid] /= norm(normals[nid])
    end

    @debug("Calculate contact pairs")
    for slave_element in slave_elements
        slave_element_nodes = get_connectivity(slave_element)
        c1, c2 = get_connectivity(slave_element)
        X1 = slave_element("geometry", time)
        if haskey(slave_element, "displacement")
            u1 = slave_element("displacement", time)
            x1 = ((Xi+ui for (Xi,ui) in zip(X1,u1))...,)
        else
            x1 = X1
        end
        n1 = ((normals[i] for i in slave_element_nodes)...,)
        master_elements = Set{Element}()
        for master_element in get_master_elements(problem)
            master_element_nodes = get_connectivity(master_element)
            cm1, cm2 = get_connectivity(master_element)
            X2 = master_element("geometry", time)
            if haskey(master_element, "displacement")
                u2 = master_element("displacement", time)
                x2 = ((Xi+ui for (Xi,ui) in zip(X2,u2))...,)
            else
                x2 = X2
            end

            x1_midpoint = 1/2*(x1[1]+x1[2])
            x2_midpoint = 1/2*(x2[1]+x2[2])
            distance = norm(x2_midpoint - x1_midpoint)
            distance > 3.0*norm(x1[2]-x1[1]) && continue

            xi1a = project_from_master_to_slave(Val{:Seg2}, x2[1], x1[1], x1[2], n1[1], n1[2])
            xi1b = project_from_master_to_slave(Val{:Seg2}, x2[2], x1[1], x1[2], n1[1], n1[2])
            xi1a = clamp(xi1a, -1.0, 1.0)
            xi1b = clamp(xi1b, -1.0, 1.0)
            l = 1/2*abs(xi1b-xi1a)
            isapprox(l, 0.0) && continue
            push!(master_elements, master_element)
        end
        n = length(master_elements)
        @debug("Contact slave element $(slave_element.id) has $n master elements.")
        props.contact_pairs[slave_element] = master_elements
    end
end

include("lasse.jl")
include("nasse.jl")
include("passe.jl")

function FEMBase.assemble_elements!(problem::Problem{Contact2DAD}, assembly::Assembly,
                                    elements::Vector{Element{Seg2}}, time::Float64)
    n = problem.properties.algorithm
    if n == 1
        assemble_interface_1!(problem, assembly, elements, time)
    elseif n == 2
        assemble_interface_2!(problem, assembly, elements, time)
    elseif n == 3
        assemble_interface_3!(problem, assembly, elements, time)
    end
    return nothing
end
