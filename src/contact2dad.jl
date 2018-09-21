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
    return Contact2DAD(master_elements, rotate_normals, dual_basis, iteration,
                       contact_state_in_first_iteration, always_in_contact,
                       use_scaling, alpha, beta, nadj_nodes, scaling_factors,
                       print_summary, contact_pairs)
end

function FEMBase.add_elements!(::Problem{Contact2DAD}, ::Any)
    error("use `add_slave_elements!` and `add_master_elements!` to add ",
          "elements to the Mortar2D problem.")
end

function FEMBase.add_slave_elements!(problem::Problem{Contact2DAD}, elements)
    for element in elements
        push!(problem.elements, element)
    end
end

function FEMBase.add_master_elements!(problem::Problem{Contact2DAD}, elements)
    for element in elements
        push!(problem.properties.master_elements, element)
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

function get_slave_dofs(problem::Problem{Contact2DAD})
    dofs = Int64[]
    for element in get_slave_elements(problem)
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

function get_slave_nodes(problem::Problem{Contact2DAD})
    nodes = Set{Int64}()
    for element in get_slave_elements(problem)
        for j in get_connectivity(element)
            push!(nodes, j)
        end
    end
    return nodes
end

function get_master_dofs(problem::Problem{Contact2DAD})
    dofs = Int64[]
    for element in get_master_elements(problem)
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

function get_master_nodes(problem::Problem{Contact2DAD})
    nodes = Set{Int64}()
    for element in get_master_elements(problem)
        for j in get_connectivity(element)
            push!(nodes, j)
        end
    end
    return nodes
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

function FEMBase.assemble_elements!(problem::Problem{Contact2DAD}, assembly::Assembly,
                                    elements::Vector{Element{Seg2}}, time::Float64)

    props = problem.properties
    props.iteration += 1

    slave_elements = get_slave_elements(problem)
    X = problem("geometry", time)
    S = get_slave_nodes(problem)
    M = get_master_nodes(problem)

    if props.iteration == 1
        assemble_elements_preprocess!(problem, assembly, elements, time)
    end

    n_evaluations = 0

    function calculate_interface(x::Vector{T}) where {T}

        dim = 2
        ndofs = Int(length(x)/2)
        nnodes = Int(ndofs/dim)
        u = Dict(j => [x[dim*(j-1)+i] for i=1:dim] for j=1:nnodes)
        la = Dict(j => [x[dim*(j-1)+i+ndofs] for i=1:dim] for j=1:nnodes)
        x = Dict(j => X[j]+u[j] for j in union(S, M))

        fc = zeros(T, 2, nnodes)
        gap = zeros(T, 2, nnodes)
        C = zeros(T, 2, nnodes)

        # 1. update nodal normals for slave elements
        normals = Dict{Int64, Vector{T}}(nid => zeros(T, 2) for nid in S)
        tangents = Dict{Int64, Vector{T}}(nid => zeros(T, 2) for nid in S)
        coeff = props.rotate_normals ? -1.0 : 1.0
        for element in slave_elements
            a, b = conn = get_connectivity(element)
            x1, x2 = x[a], x[b]
            t = t1, t2 = coeff * (x2-x1) / norm(x2-x1)
            n = [-t2, t1]
            for c in conn
                normals[c] += n
                tangents[c] += t
            end
        end
        for nid in keys(normals)
            normals[nid] /= norm(normals[nid])
            tangents[nid] /= norm(tangents[nid])
        end

        update!(slave_elements, "normal", time => normals)
        update!(slave_elements, "tangent", time => tangents)

        alpha = empty!(problem.properties.alpha)
        beta = empty!(problem.properties.beta)
        nadj_nodes = empty!(problem.properties.nadj_nodes)

        # 2. loop all slave elements
        for slave_element in slave_elements

            slave_element_nodes = get_connectivity(slave_element)
            c1, c2 = get_connectivity(slave_element)
            X1 = (X[c1], X[c2])
            u1 = (u[c1], u[c2])
            x1 = (x[c1], x[c2])
            la1 = (la[c1], la[c2])
            n1 = (normals[c1], normals[c2])
            t1 = (tangents[c1], tangents[c2])
            nnodes = size(slave_element, 2)

            Ae = Matrix(1.0I, nnodes, nnodes)

            if problem.properties.dual_basis # construct dual basis

                De = zeros(nnodes, nnodes)
                Me = zeros(nnodes, nnodes)
                ae = zeros(nnodes)
                be = zeros(nnodes)

                for ip in get_integration_points(slave_element, 3)
                    dN = get_dbasis(slave_element, ip, time)
                    j = sum([kron(dN[:,i], x1[i]') for i=1:length(x1)])
                    detj = norm(j)
                    w = ip.weight*detj
                    N1 = slave_element(ip, time)
                    ae += w*vec(N1)/detj
                end

                nsegments = 0

                for master_element in get_master_elements(problem, slave_element)

                    master_element_nodes = get_connectivity(master_element)
                    cm1, cm2 = get_connectivity(master_element)
                    X2 = (X[cm1], X[cm2])
                    u2 = (u[cm1], u[cm2])
                    x2 = (x[cm1], x[cm2])

                    # calculate segmentation: we care only about endpoints
                    xi1a = project_from_master_to_slave(Val{:Seg2}, x2[1], x1[1], x1[2], n1[1], n1[2])
                    xi1b = project_from_master_to_slave(Val{:Seg2}, x2[2], x1[1], x1[2], n1[1], n1[2])
                    xi1a = clamp(xi1a, -1.0, 1.0)
                    xi1b = clamp(xi1b, -1.0, 1.0)
                    l = 1/2*abs(xi1b-xi1a)
                    isapprox(l, 0.0) && continue # no contribution in this master element

                    nsegments += 1
                    for ip in get_integration_points(slave_element, 3)
                        # jacobian of slave element in deformed state
                        dN = get_dbasis(slave_element, ip, time)
                        j = sum([kron(dN[:,i], x1[i]') for i=1:length(x1)])
                        detj = norm(j)
                        w = ip.weight*norm(j)*l
                        xi = first(ip.coords)
                        xi_s = 1/2*(1-xi)*xi1a + 1/2*(1+xi)*xi1b
                        N1 = get_basis(slave_element, xi_s, time)
                        De += w*Matrix(Diagonal(vec(N1)))
                        Me += w*N1'*N1
                        be += w*vec(N1)/detj
                    end
                end

                for (i, j) in enumerate(get_connectivity(slave_element))
                    alpha[j] = get(alpha, j, 0.0) + ae[i]
                    beta[j] = get(beta, j, 0.0) + be[i]
                    nadj_nodes[j] = get(nadj_nodes, j, 0) + 1
                end

                if nsegments == 0
                    continue
                end

                Ae = De*inv(Me)

            end

            # 3. loop all master elements
            for master_element in get_master_elements(problem, slave_element)

                master_element_nodes = get_connectivity(master_element)
                cm1, cm2 = get_connectivity(master_element)
                X2 = (X[cm1], X[cm2])
                u2 = (u[cm1], u[cm2])
                x2 = (x[cm1], x[cm2])

                # calculate segmentation: we care only about endpoints
                xi1a = project_from_master_to_slave(Val{:Seg2}, x2[1], x1[1], x1[2], n1[1], n1[2])
                xi1b = project_from_master_to_slave(Val{:Seg2}, x2[2], x1[1], x1[2], n1[1], n1[2])
                xi1a = clamp(xi1a, -1.0, 1.0)
                xi1b = clamp(xi1b, -1.0, 1.0)
                l = 1/2*abs(xi1b-xi1a)
                isapprox(l, 0.0) && continue # no contribution in this master element

                slave_dofs = get_gdofs(problem, slave_element)
                master_dofs = get_gdofs(problem, master_element)

                # 4. loop integration points of segment
                for ip in get_integration_points(slave_element, 3)
                    # jacobian of slave element in deformed state
                    dN = get_dbasis(slave_element, ip, time)
                    j = sum([kron(dN[:,i], x1[i]') for i=1:length(x1)])
                    w = ip.weight*norm(j)*l

                    # project gauss point from slave element to master element
                    xi = first(ip.coords)
                    xi_s = 1/2*(1-xi)*xi1a + 1/2*(1+xi)*xi1b
                    N1 = vec(get_basis(slave_element, xi_s, time))
                    x_s = interpolate(N1, x1) # coordinate in gauss point
                    n_s = interpolate(N1, n1) # normal direction in gauss point
                    t_s = interpolate(N1, t1) # tangent direction in gauss point
                    #xi_m = project_from_slave_to_master_ad(master_element, x_s, n_s, x2)
                    xi_m = project_from_slave_to_master(Val{:Seg2}, x_s, n_s, x2[1], x2[2])
                    N2 = vec(get_basis(master_element, xi_m, time))
                    x_m = interpolate(N2, x2)
                    Phi = Ae*N1

                    la_s = interpolate(Phi, la1) # traction force in gauss point
                    gn = -dot(n_s, x_s - x_m) # normal gap

                    fc[:,slave_element_nodes] += w*la_s*N1'
                    fc[:,master_element_nodes] -= w*la_s*N2'
                    gap[1,slave_element_nodes] += w*gn*Phi
                    #gap[1,slave_element_nodes] += w*gn*N1'

                end # done integrating segment

            end # master elements done

        end # slave elements done

        # at this point we have calculated contact force fc and gap for all slave elements.
        # next task is to find out are they in contact or not and remove inactive nodes

        state = problem.properties.contact_state_in_first_iteration
        if problem.properties.iteration == 1
            @info("First contact iteration, initial contact state = $state")
            if state == :AUTO
                avg_gap = ForwardDiff.value(mean([gap[1, j] for j in S]))
                std_gap = ForwardDiff.value(std([gap[1, j] for j in S]))
                if (avg_gap < 1.0e-3) && (std_gap < 1.0e-6)
                    state = :ACTIVE
                else
                    state = :UNKNOWN
                end
                @info("Average weighted gap = $avg_gap, std gap = $std_gap, automatically determined contact state = $state")
            end
        end

        is_active = Dict{Int, Bool}()
        condition = Dict()

        for j in S
            if j in props.always_in_contact
                is_active[j] = true
                continue
            end
            lan = dot(normals[j], la[j])
            condition[j] = ForwardDiff.value(lan - gap[1, j])
            is_active[j] = condition[j] > 0
        end

        if problem.properties.iteration == 1 && state == :ACTIVE
            for j in S
                is_active[j] = true
                condition[j] = 0.0
            end
        end

        if problem.properties.iteration == 1 && state == :INACTIVE
            for j in S
                is_active[j] = false
                condition[j] = 0.0
            end
        end

        if problem.properties.print_summary
            @info("Summary of nodes")
            if problem.properties.iteration == 1 && state == :ACTIVE
                @info("Contact state in first iteration determined to be ACTIVE")
            elseif problem.properties.iteration == 1 && state == :INACTIVE
                @info("Contact state in first iteration determined to be INACTIVE")
            else
                for j in sort(collect(keys(is_active)))
                    n = map(ForwardDiff.value, normals[j])
                    @info("$j, c=$(condition[j]), s=$(is_active[j]), n=$n")
                end
            end
        end

        for j in S

            if is_active[j]
                lan = dot(normals[j], la[j])
                lat = dot(tangents[j], la[j])
                C[1,j] -= gap[1, j]
                C[2,j] += lat
            else
                C[:,j] = la[j]
            end

        end

        # Apply scaling to contact virtual work and contact constraints
        if problem.properties.use_scaling
            scaling = empty!(problem.properties.scaling_factors)
            for j in S
                isapprox(beta[j], 0.0) && continue
                scaling[j] = alpha[j] / (nadj_nodes[j] * beta[j])
                C[1,j] *= scaling[j]
                fc[:,j] *= scaling[j]
            end
            update!(slave_elements, "contact scaling", time => scaling)
        end

        n_evaluations += 1

        return vec([fc C])

    end

    # Next, construct vector x = [u; la] what we are going to use as input
    S = get_slave_nodes(problem)
    M = get_master_nodes(problem)
    max_node_id = max(maximum(S), maximum(M))
    n_slave_nodes = length(S)
    n_master_nodes = length(M)
    @info("Interface statistics", max_node_id, n_slave_nodes, n_master_nodes)

    ndofs = max_node_id * 2
    if max_node_id > n_slave_nodes+n_master_nodes
        warn("highly unefficient differentiation coming.")
    end

    if length(problem.assembly.u) == 0
        error("2d autodiff contact problem: initialize problem.assembly.u & la before solution")
    end
    x = [problem.assembly.u[1:ndofs]; problem.assembly.la[1:ndofs]]
    chunk_size = length(x)
    @info("ndofs = $ndofs, chunk size = $chunk_size")
    chunk = ForwardDiff.Chunk{chunk_size}()
    cfg = ForwardDiff.JacobianConfig(calculate_interface, x, chunk)

    result = DiffResults.JacobianResult(x)
    result = ForwardDiff.jacobian!(result, calculate_interface, x, cfg)
    #result = ForwardDiff.jacobian!(result, calculate_interface, x)
    A = DiffResults.jacobian(result)
    b = DiffResults.value(result)
    @info("Interface was evaluated $n_evaluations times during automatic differentiation.")

    # x doesn't mean deformed configuration here
    # x = [problem.assembly.u; problem.assembly.la]
    # if length(x) == 0
    #     error("2d autodiff contact problem: initialize problem.assembly.u & la before solution")
    # end
    # A = ForwardDiff.jacobian(calculate_interface, x)
    # b = calculate_interface(x)

    A = sparse(A)
    b = -sparse(b)
    SparseArrays.droptol!(A, 1.0e-9)
    SparseArrays.droptol!(b, 1.0e-9)

    ndofs = round(Int, length(x)/2)
    K = A[1:ndofs,1:ndofs]
    C1 = A[1:ndofs,ndofs+1:end]
    C2 = A[ndofs+1:end,1:ndofs]
    D = A[ndofs+1:end,ndofs+1:end]
    f = b[1:ndofs]
    g = b[ndofs+1:end]

    f += C1*problem.assembly.la[1:ndofs]
    g += D*problem.assembly.la[1:ndofs]

    # if !haskey(problem, "contact force")
    #     problem.fields["contact force"] = field(time => f)
    # else
    #     update!(problem.fields["contact force"], time => f)
    # end
    # fc = problem.fields["contact force"]
    # if length(fc) > 1
    #     # kick in generalized alpha rule for time integration
    #     alpha = 0.5
    #     @info("Applying Generalized alpha time integration")
    #     K = (1-alpha)*K
    #     C1 = (1-alpha)*C1
    #     f = alpha*fc[end-1].data
    # end

    problem.assembly.K = K
    problem.assembly.C1 = copy(transpose(C1))
    problem.assembly.C2 = C2
    problem.assembly.D = D
    problem.assembly.f = sparse(f)
    problem.assembly.g = sparse(g)

end
