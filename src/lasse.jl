# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

function assemble_interface_1!(problem::Problem{Contact2DAD}, assembly::Assembly,
                               elements::Vector{Element{Seg2}}, time::Float64)

    props = problem.properties
    props.iteration += 1
    if props.iteration == 1
        assemble_elements_preprocess!(problem, assembly, elements, time)
    end

    n_evaluations = 0

    function F!(y::Vector{T}, x::Vector{T}) where {T}

        slave_elements = get_slave_elements(problem)
        X = problem("geometry", time)
        S = get_slave_nodes(problem)
        M = get_master_nodes(problem)

        u = Dict(j => [x[2*(j-1)+i] for i=1:2] for j in union(S, M))
        la = Dict(j => [x[2*(j-1)+i+ndofs] for i=1:2] for j in S)
        x = Dict(j => X[j]+u[j] for j in union(S, M))

        f = Dict(j => zeros(T, 2) for j in union(S, M))
        g = Dict(j => zeros(T, 2) for j in S)
        n = Dict(j => zeros(T, 2) for j in S)
        t = Dict(j => zeros(T, 2) for j in S)

        gn = Dict(j => zero(T) for j in S)
        lan = Dict(j => zero(T) for j in S)
        lat = Dict(j => zero(T) for j in S)

        # 1. update nodal normals for slave elements
        coeff = props.rotate_normals ? -1.0 : 1.0
        for element in slave_elements
            a, b = conn = get_connectivity(element)
            x1, x2 = x[a], x[b]
            t1, t2 = coeff * (x2-x1) / norm(x2-x1)
            n1, n2 = -t2, t1
            for c in conn
                n[c] += [n1, n2]
                t[c] += [t1, t2]
            end
        end
        for j in S
            n[j] /= norm(n[j])
            t[j] /= norm(t[j])
        end

        # 2. loop all slave elements
        for slave_element in slave_elements

            cs1, cs2 = get_connectivity(slave_element)

            Ae = zeros(T, 2, 2)
            if problem.properties.dual_basis # construct dual basis
                calculate_dual_basis_coefficients!(Ae, problem, slave_element, x, n)
            else
                Ae[:,:] .= Matrix(1.0I, 2, 2)
            end

            # 3. loop all master elements
            for master_element in get_master_elements(problem, slave_element)

                cm1, cm2 = get_connectivity(master_element)

                # calculate segmentation: we care only about endpoints
                xi1a = project_from_master_to_slave(Val{:Seg2}, x[cm1], x[cs1], x[cs2], n[cs1], n[cs2])
                xi1b = project_from_master_to_slave(Val{:Seg2}, x[cm2], x[cs1], x[cs2], n[cs1], n[cs2])
                xi1a = clamp(xi1a, -1.0, 1.0)
                xi1b = clamp(xi1b, -1.0, 1.0)
                l = 1/2*abs(xi1b-xi1a)
                isapprox(l, 0.0) && continue # no contribution in this master element

                # 4. loop integration points of segment
                for (weight, xi) in quadrature_points
                    # jacobian of slave element in deformed state
                    detj = norm(0.5*(x[cs2]-x[cs1]))
                    w = weight*detj*l

                    # project gauss point from slave element to master element
                    xi_s = 0.5*(1.0-xi)*xi1a + 0.5*(1.0+xi)*xi1b
                    N1 = (0.5*(1.0-xi_s), 0.5*(1.0+xi_s))
                    x_s = interpolate(N1, (x[cs1], x[cs2])) # coordinate in gauss point
                    n_s = interpolate(N1, (n[cs1], n[cs2])) # normal direction in gauss point
                    t_s = interpolate(N1, (t[cs1], t[cs2])) # tangent direction in gauss point
                    normalize = true
                    if normalize
                        n_s ./= norm(n_s)
                        t_s ./= norm(t_s)
                    end
                    Phi = Ae*collect(N1)
                    la_s = interpolate(Phi, (la[cs1], la[cs2])) # traction force in gauss point

                    xi_m = project_from_slave_to_master(Val{:Seg2}, x_s, n_s, x[cm1], x[cm2])
                    N2 = (0.5*(1.0-xi_m), 0.5*(1.0+xi_m))
                    x_m = interpolate(N2, (x[cm1], x[cm2]))

                    slave_element_nodes = get_connectivity(slave_element)
                    master_element_nodes = get_connectivity(master_element)
                    for (i,j) in enumerate(slave_element_nodes)
                        f[j] += w*la_s*N1[i]
                    end
                    for (i,j) in enumerate(master_element_nodes)
                        f[j] -= w*la_s*N2[i]
                    end
                    for (i,j) in enumerate(slave_element_nodes)
                        g[j] += w*(x_m-x_s)*Phi[i]
                        gn[j] += w*dot(n_s, x_m-x_s)*Phi[i]
                        lan[j] += w*dot(n_s, la_s)*Phi[i]
                        lat[j] += w*dot(n_s, la_s)*Phi[i]
                    end

                end # done integrating segment

            end # master elements done

        end

        # at this point we have calculated contact force fc and gap for all slave elements.
        # next task is to find out are they in contact or not and remove inactive nodes

        is_active = Dict{Int64, Bool}()
        complementarity_condition = Dict{Int64, Float64}()
        # lan = Dict(j => dot(n[j], la[j]) for j in S)
        # lat = Dict(j => dot(t[j], la[j]) for j in S)
        # gn = Dict(j => dot(n[j], g[:,j]) for j in S)

        state = problem.properties.contact_state_in_first_iteration
        if problem.properties.iteration == 1
            @info("First contact iteration, initial contact state = $state")
            if state == :AUTO
                avg_gap = value(mean(gn[j] for j in S))
                std_gap = value(std(gn[j] for j in S))
                if (avg_gap < 1.0e-3) && (std_gap < 1.0e-6)
                    state = :ACTIVE
                else
                    state = :UNKNOWN
                end
                @info("Average weighted gap = $avg_gap, std gap = $std_gap, automatically determined contact state = $state")
            end
        end

        for j in S
            if j in props.always_in_contact
                @info("Set node $j to be always in contact.")
                is_active[j] = true
                complementarity_condition[j] = 1.0
                continue
            end
            complementarity_condition[j] = value(lan[j] - gn[j])
            is_active[j] = complementarity_condition[j] > 0
        end

        if problem.properties.iteration == 1
            for j in S
                if state == :ACTIVE
                    is_active[j] = true
                    complementarity_condition[j] = 1.0
                elseif state == :INACTIVE
                    is_active[j] = false
                    complementarity_condition[j] = -1.0
                end
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
                    Cn = round(complementarity_condition[j], digits=3)
                    normal = round.(map(value, n[j]), digits=3)
                    lambda_normal = round.(map(value, lan[j]), digits=3)
                    gap_normal = round.(map(value, gn[j]), digits=3)
                    status = is_active[j] ? "ACTIVE" : "INACTIVE"
                    @info("n($j)=$normal, " *
                          "(g⋅n)($j)=$gap_normal, " *
                          "(λ⋅n)($j)=$lambda_normal, " *
                          "(C⋅n)($j)=$Cn, status = $status")
                end
            end
        end

        if problem.properties.use_scaling

            # Apply scaling to contact virtual work and contact constraints

            alpha = Dict(j => zero(T) for j in S)
            beta = Dict(j => zero(T) for j in S)
            scaling = Dict(j => zero(T) for j in S)
            nadj_nodes = Dict(j => zero(Int64) for j in S)

            for slave_element in slave_elements
                cs1, cs2 = get_connectivity(slave_element)
                ae = zeros(T, 2)
                be = zeros(T, 2)
                for (weight, xi) in quadrature_points
                    detj = norm(0.5*(x[cs2]-x[cs1]))
                    w = weight*detj
                    N1 = [0.5*(1.0-xi), 0.5*(1.0+xi)]
                    ae += w*vec(N1)/detj
                end
                for master_element in get_master_elements(problem, slave_element)
                    cm1, cm2 = get_connectivity(master_element)
                    xi1a = project_from_master_to_slave(Val{:Seg2}, x[cm1], x[cs1], x[cs2], n[cs1], n[cs2])
                    xi1b = project_from_master_to_slave(Val{:Seg2}, x[cm2], x[cs1], x[cs2], n[cs1], n[cs2])
                    xi1a = clamp(xi1a, -1.0, 1.0)
                    xi1b = clamp(xi1b, -1.0, 1.0)
                    l = 1/2*abs(xi1b-xi1a)
                    isapprox(l, 0.0) && continue
                    for (weight, xi) in quadrature_points
                        detj = norm(0.5*(x[cs2]-x[cs1]))
                        w = weight*detj*l
                        N1 = [0.5*(1.0-xi), 0.5*(1.0+xi)]
                        be += w*N1/detj
                    end
                end
                for (i, j) in enumerate(get_connectivity(slave_element))
                    alpha[j] += ae[i]
                    beta[j] += be[i]
                    nadj_nodes[j] += 1
                end
            end

            for j in S
                isapprox(beta[j], 0.0) && continue
                scaling[j] = alpha[j] / (nadj_nodes[j] * beta[j])
                gn[j] *= scaling[j]
                f[j] *= scaling[j]
            end

            update!(slave_elements, "contact scaling", time => scaling)

        end

        # Contact force
        for j in union(S, M)
            dofs = [2*(j-1)+i for i=1:2]
            y[dofs] .= f[j]
        end

        # Contact constraints
        for j in S
            offset = 2*(length(S)+length(M))
            dofs = dof1, dof2 = [2*(j-1)+offset+i for i=1:2]
            if is_active[j]
                y[dof1] = -gn[j]
                y[dof2] = lat[j]
            else
                y[dofs] = la[j]
            end
        end

        n_evaluations += 1

        return nothing

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
        @warn("Highly unefficient differentiation coming. To make differentiation " *
              "more efficient, reorder nodes so that interface nodes comes first")
    end

    if length(problem.assembly.u) == 0
        error("2d autodiff contact problem: initialize problem.assembly.u & la before solution")
    end

    x = [problem.assembly.u[1:ndofs]; problem.assembly.la[1:ndofs]]
    y = zeros(length(x))
    chunk_size = length(x)
    @info("ndofs = $ndofs, chunk size = $chunk_size")
    chunk = ForwardDiff.Chunk{chunk_size}()
    cfg = ForwardDiff.JacobianConfig(F!, y, x, chunk)
    result = DiffResults.JacobianResult(y, x)
    ForwardDiff.jacobian!(result, F!, y, x, cfg)
    A = sparse(DiffResults.jacobian(result))
    b = -sparse(DiffResults.value(result))
    SparseArrays.droptol!(A, 1.0e-9)
    SparseArrays.droptol!(b, 1.0e-9)
    @info("Interface was evaluated $n_evaluations times during automatic differentiation.")

    ndofs = round(Int, length(x)/2)
    K = A[1:ndofs,1:ndofs]
    C1 = A[1:ndofs,ndofs+1:end]
    C2 = A[ndofs+1:end,1:ndofs]
    D = A[ndofs+1:end,ndofs+1:end]
    f = b[1:ndofs] + C1*problem.assembly.la[1:ndofs]
    g = b[ndofs+1:end] + D*problem.assembly.la[1:ndofs]

    problem.assembly.K = K
    problem.assembly.C1 = copy(transpose(C1))
    problem.assembly.C2 = C2
    problem.assembly.D = D
    problem.assembly.f = sparse(f)
    problem.assembly.g = sparse(g)

    return nothing

end
