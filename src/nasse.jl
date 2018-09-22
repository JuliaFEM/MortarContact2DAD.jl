# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2DAD.jl/blob/master/LICENSE

quadrature_points = FEMBase.get_quadrature_points(Val{:GLSEG5})

function nodal_assembly_factory(node_id, problem, time)

    X = problem("geometry", time)
    slave_elements = filter(element -> in(node_id, get_connectivity(element)), get_slave_elements(problem))
    master_elements = union([get_master_elements(problem, slave_element) for slave_element in slave_elements]...)
    slave_nodes = collect(get_nodes(slave_elements))
    S = slave_nodes = [node_id; setdiff(slave_nodes, node_id)]
    M = master_nodes = collect(get_nodes(master_elements))
    slave_dofs = [2*(j-1)+i for j in slave_nodes for i=1:2]
    master_dofs = [2*(j-1)+i for j in master_nodes for i=1:2]
    n_slave_elements = length(slave_elements)
    n_master_elements = length(slave_elements)
    n_slave_nodes = length(slave_nodes)
    n_master_nodes = length(master_nodes)
    n_slave_dofs = length(slave_dofs)
    n_master_dofs = length(master_dofs)
    n_dofs = n_slave_dofs+n_master_dofs
    @info("Contact interface statistics:")
    @info("$n_slave_elements slave elements connecting to node $node_id.")
    @info("$n_master_elements master elements pairing with slave elements.")
    @info("$n_slave_nodes slave nodes: $slave_nodes")
    @info("$n_slave_dofs slave dofs: $slave_dofs")
    @info("$n_master_nodes master nodes: $master_nodes")
    @info("$n_master_dofs master_dofs: $master_dofs")
    @info("slave node id in calculation: $node_id")
    nodemap = Dict(j => i for (i, j) in enumerate(union(slave_nodes, master_nodes)))
    @info("nodes mapping inside algorithm: $nodemap")

    function F!(y::AbstractArray, x::Vector{T}) where {T}

        # unpack input data to dictionaries for easier usage
        # the order is: first slave nodes (with node_id first)
        # then master nodes, and lastly Lagrange multiplier of slave nodes

        u = Dict(j => [x[2*(nodemap[j]-1)+i] for i=1:2] for j in union(S, M))
        la = Dict(j => [x[2*(nodemap[j]-1)+i+n_dofs] for i=1:2] for j in S)
        x = Dict(j => X[j]+u[j] for j in union(S, M))
        n = Dict(j => zeros(T, 2) for j in S)
        t = Dict(j => zeros(T, 2) for j in S)
        g = Dict(j => zeros(T, 2) for j in S)
        f = Dict(j => zeros(T, 2) for j in union(S,M))

        # 1. update nodal normals for slave elements
        coeff = problem.properties.rotate_normals ? -1.0 : 1.0
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

        alpha = empty!(problem.properties.alpha)
        beta = empty!(problem.properties.beta)
        nadj_nodes = empty!(problem.properties.nadj_nodes)

        # 2. loop all slave elements
        for slave_element in slave_elements

            cs1, cs2 = get_connectivity(slave_element)

            if problem.properties.dual_basis # construct dual basis

                De = zeros(T, 2, 2)
                Me = zeros(T, 2, 2)
                ae = zeros(T, 2)
                be = zeros(T, 2)

                for (weight, xi) in quadrature_points
                    detj = norm(0.5*(x[cs2]-x[cs1]))
                    w = weight*detj
                    N1 = [0.5*(1.0-xi), 0.5*(1.0+xi)]
                    ae += w*N1/detj
                end

                nsegments = 0

                for master_element in master_elements

                    cm1, cm2 = get_connectivity(master_element)

                    # calculate segmentation
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
                        Me += w*kron(N1', N1)
                        be += w*N1/detj
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
            else
                Ae = Matrix(1.0I, 2, 2)
            end

            # Dual basis constructed, create interface

            for master_element in master_elements

                cm1, cm2 = get_connectivity(master_element)

                # calculate segmentation
                xi1a = project_from_master_to_slave(Val{:Seg2}, x[cm1], x[cs1], x[cs2], n[cs1], n[cs2])
                xi1b = project_from_master_to_slave(Val{:Seg2}, x[cm2], x[cs1], x[cs2], n[cs1], n[cs2])
                xi1a = clamp(xi1a, -1.0, 1.0)
                xi1b = clamp(xi1b, -1.0, 1.0)
                l = 1/2*abs(xi1b-xi1a)
                isapprox(l, 0.0) && continue # no contribution in this master element

                # 4. loop integration points of segment
                for (weight, xi) in quadrature_points
                    detj = norm(0.5*(x[cs2]-x[cs1]))
                    w = weight*detj*l

                    # project gauss point from slave element to master element
                    xi_s = 0.5*(1.0-xi)*xi1a + 0.5*(1.0+xi)*xi1b
                    N1 = (0.5*(1.0-xi_s), 0.5*(1.0+xi_s))
                    x1 = (x[cs1], x[cs2])
                    n1 = (n[cs1], n[cs2])
                    t1 = (t[cs1], t[cs2])
                    la1 = (la[cs1], la[cs2])
                    x_s = interpolate(N1, x1) # coordinate in gauss point
                    n_s = interpolate(N1, n1) # normal direction in gauss point
                    t_s = interpolate(N1, t1) # tangent direction in gauss point
                    Phi = Ae*collect(N1)
                    la_s = interpolate(Phi, la1) # traction force in gauss point

                    xi_m = project_from_slave_to_master(Val{:Seg2}, x_s, n_s, x[cm1], x[cm2])
                    N2 = (0.5*(1.0-xi_m), 0.5*(1.0+xi_m))
                    x2 = (x[cm1], x[cm2])
                    x_m = interpolate(N2, x2)

                    # Write contact force of slave element
                    for (i,j) in enumerate(get_connectivity(slave_element))
                        f[j] += w*la_s*N1[i]
                    end

                    # Write contact force of master element
                    for (i,j) in enumerate(get_connectivity(master_element))
                        f[j] -= w*la_s*N2[i]
                    end

                    # Write weighted gap
                    for (i,j) in enumerate(get_connectivity(slave_element))
                        g[j] += w*(x_m-x_s)*Phi[i]
                    end

                end # done integrating segment

            end # master elements done

        end # slave elements done

        # at this point we have calculated contact force fc and gap for all slave elements.
        # next task is to find out are they in contact or not and remove inactive nodes

        is_active = Dict{Int64, Bool}()
        complementarity_condition = Dict{Int64, Float64}()
        lan = Dict(j => dot(n[j], la[j]) for j in S)
        lat = Dict(j => dot(t[j], la[j]) for j in S)
        gn = Dict(j => dot(n[j], g[j]) for j in S)

        for j in S
            if j in problem.properties.always_in_contact
                is_active[j] = true
                complementarity_condition[j] = 0.0
            else
                complementarity_condition[j] = value(lan[j] - gn[j])
                is_active[j] = complementarity_condition[j] >= 0
            end
        end

        if problem.properties.iteration == 1
            state = problem.properties.contact_state_in_first_iteration
            if state == :ACTIVE
                for j in S
                    is_active[j] = true
                    complementarity_condition[j] = 1.0
                end
            elseif state == :INACTIVE
                for j in S
                    is_active[j] = false
                    complementarity_condition[j] = -1.0
                end
            end
        end

        if problem.properties.print_summary
            @info("Summary for node $node_id")
            if problem.properties.iteration == 1 && state == :ACTIVE
                @info("Contact state in first iteration determined to be ACTIVE")
            elseif problem.properties.iteration == 1 && state == :INACTIVE
                @info("Contact state in first iteration determined to be INACTIVE")
            else
                Cn = round(complementarity_condition[node_id], digits=3)
                normal = round.(map(value, n[node_id]), digits=3)
                lambda_normal = round.(map(value, lan[node_id]), digits=3)
                gap_normal = round.(map(value, gn[node_id]), digits=3)
                status = is_active[node_id] ? "ACTIVE" : "INACTIVE"
                @info("n($node_id)=$normal, " *
                      "(g⋅n)($node_id)=$gap_normal, " *
                      "(λ⋅n)($node_id)=$lambda_normal, " *
                      "(C⋅n)($node_id)=$Cn, status = $status")
            end
        end

        # First update contact virtual work part
        for j in union(slave_nodes, master_nodes)
            dofs = [2*(nodemap[j]-1)+i for i=1:2]
            @debug("write f for node $j to dofs $dofs")
            y[dofs] = f[j]
        end

        # Then update contact constraint part
        for j in slave_nodes
            dofs = dof1, dof2 = [2*(nodemap[j]-1)+n_dofs+i for i=1:2]
            @debug("write contact constraint cn($j) to dofs $dofs")
            if is_active[j]
                y[dof1] = -gn[node_id]
                y[dof2] = lat[node_id]
            else
                y[dofs] = la[node_id]
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

        return y

    end # F!

    return F!, slave_dofs, master_dofs

end


function assemble_interface_2!(problem::Problem{Contact2DAD}, assembly::Assembly,
                               elements::Vector{Element{Seg2}}, time::Float64)

    props = problem.properties
    props.iteration += 1
    if props.iteration == 1
        assemble_elements_preprocess!(problem, assembly, elements, time)
    end

    S = get_slave_nodes(problem)
    n_slave_nodes = length(S)
    @debug("Contact interface has $n_slave_nodes nodes")

    for j in sort(collect(S))
        @debug("Finding solution for node $j")
        F!, slave_dofs, master_dofs = nodal_assembly_factory(j, problem, time)
        u = problem.assembly.u
        la = problem.assembly.la
        us = u[slave_dofs]
        um = u[master_dofs]
        las = la[slave_dofs]
        x = [us; um; las]
        y = zeros(length(x))
        use_cache = true
        if use_cache
            @info("Using DiffResults to store low-level results.")
            chunk_size = length(x)
            chunk = ForwardDiff.Chunk{chunk_size}()
            cfg = ForwardDiff.JacobianConfig(F!, y, x, chunk)
            result = DiffResults.JacobianResult(y, x)
            @info("Evaluating contact interface, chunk size $chunk_size")
            dt = @elapsed result = ForwardDiff.jacobian!(result, F!, y, x, cfg)
            @info("Contact interface evaluated in $dt seconds.")
            A = DiffResults.jacobian(result)
            b = -DiffResults.value(result)
        else
            info("Not using DiffResults.")
            A = ForwardDiff.jacobian(F!, y, x)
            b = -F!(y, x)
        end
        all_dofs = union(slave_dofs, master_dofs)
        ndofs = length(all_dofs)
        @info("displacement dofs in assembly", slave_dofs', master_dofs', all_dofs')
        K = A[1:ndofs, 1:ndofs]
        C1 = A[1:ndofs, ndofs+1:end]
        C2 = A[ndofs+1:end, 1:ndofs]
        D = A[ndofs+1:end, ndofs+1:end]
        f = b[1:ndofs] + C1*la[slave_dofs]
        g = b[ndofs+1:end] + D*la[slave_dofs]

        function empty_dofs!(dofs)
            @debug("Emptying dofs $dofs from K")
            K[dofs, :] .= 0.0
            K[:, dofs] .= 0.0
            @debug("Emptying dofs $dofs from C1")
            C1[:, dofs] .= 0.0
            @debug("Emptying dofs $dofs from C2")
            C2[dofs, :] .= 0.0
            @debug("Emptying dofs $dofs from D")
            D[dofs, :] .= 0.0
            D[:, dofs] .= 0.0
            @debug("Emptying dofs $dofs from f")
            f[dofs] .= 0.0
            @debug("Emptying dofs $dofs from g")
            g[dofs] .= 0.0
            return nothing
        end

        empty_dofs!(3:length(slave_dofs))

        droptol = true
        if droptol
            K[abs.(K) .< 1.0e-9] .= 0.0
            C1[abs.(C1) .< 1.0e-9] .= 0.0
            C2[abs.(C2) .< 1.0e-9] .= 0.0
            D[abs.(D) .< 1.0e-9] .= 0.0
            f[abs.(f) .< 1.0e-9] .= 0.0
            g[abs.(g) .< 1.0e-9] .= 0.0
        end

        @info("Contact stiffnes K = ∂fc/∂u", K)
        @info("Discrete interface operator B' = ∂fc/∂λ = [0 D -M]'", C1')
        @info("Contact constraints", C2)
        @info("Tangential constraints", D)
        @info("Contact force fc", f)
        @info("Gap vector gn", g)

        @info("Assembling to K")
        add!(assembly.K, all_dofs, all_dofs, K)
        @info("Assembling to C1")
        add!(assembly.C1, slave_dofs, all_dofs, copy(transpose(C1)))
        @info("Assembling to C2")
        add!(assembly.C2, slave_dofs, all_dofs, C2)
        @info("Assembling to D")
        add!(assembly.D, slave_dofs, slave_dofs, D)
        @info("Assembling to f")
        add!(assembly.g, slave_dofs, g)
        @info("Assembling to g")
        add!(assembly.f, all_dofs, f)
    end

    return nothing

end
