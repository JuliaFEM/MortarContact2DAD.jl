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
end

function Contact2DAD()
    return Contact2DAD([], false, true, 0, :AUTO, Set(),
                        true, Dict(), Dict(), Dict(), Dict(), false)
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

function get_slave_dofs(problem::Problem{Contact2DAD})
    dofs = Int64[]
    for element in get_slave_elements(problem)
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

function get_master_dofs(problem::Problem{Contact2DAD})
    dofs = Int64[]
    for element in get_master_elements(problem)
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

""" Find segment from slave element corresponding to master element nodes.

Parameters
----------
x1_, n1_
    slave element geometry and normal direction
x2
    master element node to project onto slave

Returns
-------
xi
    dimensionless coordinate on slave corresponding to
    projected master

"""
function project_from_master_to_slave_ad(
    slave_element::Element{E}, x1_::DVTI, n1_::DVTI, x2::Vector;
    tol=1.0e-10, max_iterations=20, debug=false) where E<:MortarElements2D

    """ Multiply basis / dbasis at `xi` with field. """
    function mul(func, xi, field)
        B = func(slave_element, [xi], time)
        return sum(B[i]*field[i] for i=1:length(B))
    end

    x1(xi1) = mul(get_basis, xi1, x1_)
    dx1(xi1) = mul(get_dbasis, xi1, x1_)
    n1(xi1) = mul(get_basis, xi1, n1_)
    dn1(xi1) = mul(get_dbasis, xi1, n1_)
    cross2(a, b) = cross([a; 0], [b; 0])[3]
    R(xi1) = cross2(x1(xi1)-x2, n1(xi1))
    dR(xi1) = cross2(dx1(xi1), n1(xi1)) + cross2(x1(xi1)-x2, dn1(xi1))

    xi1 = 0.0
    xi1_next = 0.0
    dxi1 = 0.0
    for i=1:max_iterations
        dxi1 = -R(xi1)/dR(xi1)
        dxi1 = clamp.(dxi1, -0.3, 0.3)
        xi1_next = clamp.(xi1 + dxi1, -1.0, 1.0)
        if norm(xi1_next - xi1) < tol
            return xi1_next
        end
        if debug
            @info("Debug data for iteration $i")
            @info("xi1 = $(value(xi1))")
            @info("R(xi1) = $(value(R(xi1)))")
            @info("dR(xi1) = $(value(dR(xi1)))")
            @info("dxi1 = $(value(dxi1))")
            @info("norm = $(value(norm(xi1_next - xi1)))")
            @info("xi1_next = $(value(xi1_next))")
        end
        xi1 = xi1_next
    end

    @info("Projection algorithm failed with the following data:")
    a, b = x1_.data
    @info("x11 = $(map(value, a))")
    @info("x12 = $(map(value, b))")
    a, b = n1_.data
    @info("n11 = $(map(value, a))")
    @info("n12 = $(map(value, b))")
    @info("x2 = $(map(value, x2))")
    @info("xi1 = $(value(xi1)), dxi1 = $(value(dxi1))")
    @info("-R(xi1) = $(value(-R(xi1)))")
    @info("dR(xi1) = $(value(dR(xi1)))")
    error("find projection from master to slave: did not converge")

end

function project_from_slave_to_master_ad(
    master_element::Element{E}, x1, n1, x2_;
    tol=1.0e-10, max_iterations=20) where E<:MortarElements2D

    x2(xi2) = interpolate(vec(get_basis(master_element, [xi2], time)), x2_)
    dx2(xi2) = interpolate(vec(get_dbasis(master_element, [xi2], time)), x2_)
    cross2(a, b) = cross([a; 0], [b; 0])[3]
    R(xi2) = cross2(x2(xi2)-x1, n1)
    dR(xi2) = cross2(dx2(xi2), n1)

    xi2 = 0.0
    dxi2 = 0.0
    for i=1:max_iterations
        dxi2 = -R(xi2) / dR(xi2)
        xi2 += dxi2
        if norm(dxi2) < tol
            return xi2
        end
    end

    error("find projection from slave to master: did not converge, last val: $xi2 and $dxi2")

end


"""
Frictionless 2d finite sliding contact with forwarddiff.

true/false flags: finite_sliding, friction, use_forwarddiff
"""
function FEMBase.assemble_elements!(problem::Problem{Contact2DAD}, assembly::Assembly,
                                    elements::Vector{Element{Seg2}}, time::Float64)

    props = problem.properties
    props.iteration += 1

    field_dim = get_unknown_field_dimension(problem)
    field_name = get_parent_field_name(problem)
    slave_elements = get_slave_elements(problem)

    function calculate_interface(x::Vector)

        ndofs = round(Int, length(x)/2)
        nnodes = round(Int, ndofs/field_dim)
        u = reshape(x[1:ndofs], field_dim, nnodes)
        la = reshape(x[ndofs+1:end], field_dim, nnodes)
        fc = zero(u)
        gap = zero(u)
        C = zero(la)
        S = Set{Int64}()

        # 1. update nodal normals for slave elements
        Q = [0.0 -1.0; 1.0 0.0]
        normals = zero(u)
        for element in slave_elements
            conn = get_connectivity(element)
            push!(S, conn...)
            gdofs = get_gdofs(problem, element)
            X_el = element("geometry", time)
            x_el = tuple( (X_el[i] + u[:,j] for (i,j) in enumerate(conn))... )
            #=
            for ip in get_integration_points(element, 3)
                dN = get_dbasis(element, ip, time)
                N = element(ip, time)
                t = sum([kron(dN[:,i], x_el[i]') for i=1:length(x_el)])
                normals[:, conn] += ip.weight*Q*t'*N
            end
            =#
            dN = get_dbasis(element, [0.0], time)
            t = sum([kron(dN[:,i], x_el[i]') for i=1:length(x_el)])
            n = Q*t'
            n /= norm(n)
            for c in conn
                normals[:,c] += n
            end
        end
        for i in 1:size(normals,2)
            normals[:,i] /= norm(normals[:,i])
        end
        # swap element normals in 2d if they point to inside of body
        if props.rotate_normals
            for i=1:size(normals,2)
                normals[:,i] = -normals[:,i]
            end
        end
        normals2 = Dict()
        for j in S
            normals2[j] = normals[:,j]
        end
        update!(slave_elements, "normal", time => normals2)

        alpha = empty!(problem.properties.alpha)
        beta = empty!(problem.properties.beta)
        nadj_nodes = empty!(problem.properties.nadj_nodes)

        # 2. loop all slave elements
        for slave_element in slave_elements

            slave_element_nodes = get_connectivity(slave_element)
            X1 = slave_element("geometry", time)
            u1 = ((u[:,i] for i in slave_element_nodes)...,)
            x1 = ((Xi+ui for (Xi,ui) in zip(X1,u1))...,)
            la1 = ((la[:,i] for i in slave_element_nodes)...,)
            n1 = ((normals[:,i] for i in slave_element_nodes)...,)
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

                for master_element in get_master_elements(problem)

                    master_element_nodes = get_connectivity(master_element)
                    X2 = master_element("geometry", time)
                    u2 = ((u[:,i] for i in master_element_nodes)...,)
                    x2 = ((Xi+ui for (Xi,ui) in zip(X2,u2))...,)

                    # calculate segmentation: we care only about endpoints
                    xi1a = project_from_master_to_slave_ad(slave_element, field(x1), field(n1), x2[1])
                    xi1b = project_from_master_to_slave_ad(slave_element, field(x1), field(n1), x2[2])
                    xi1 = clamp.([xi1a; xi1b], -1.0, 1.0)
                    l = 1/2*abs(xi1[2]-xi1[1])
                    isapprox(l, 0.0) && continue # no contribution in this master element

                    nsegments += 1
                    for ip in get_integration_points(slave_element, 3)
                        # jacobian of slave element in deformed state
                        dN = get_dbasis(slave_element, ip, time)
                        j = sum([kron(dN[:,i], x1[i]') for i=1:length(x1)])
                        detj = norm(j)
                        w = ip.weight*norm(j)*l
                        xi = ip.coords[1]
                        xi_s = dot([1/2*(1-xi); 1/2*(1+xi)], xi1)
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
            for master_element in get_master_elements(problem)

                master_element_nodes = get_connectivity(master_element)
                X2 = master_element("geometry", time)
                u2 = ((u[:,i] for i in master_element_nodes)...,)
                x2 = ((Xi+ui for (Xi,ui) in zip(X2,u2))...,)

                #x1_midpoint = 1/2*(x1[1]+x1[2])
                #x2_midpoint = 1/2*(x2[1]+x2[2])
                #distance = ForwardDiff.get_value(norm(x2_midpoint - x1_midpoint))
                #distance > props.maximum_distance && continue

                # calculate segmentation: we care only about endpoints
                xi1a = project_from_master_to_slave_ad(slave_element, field(x1), field(n1), x2[1])
                xi1b = project_from_master_to_slave_ad(slave_element, field(x1), field(n1), x2[2])
                xi1 = clamp.([xi1a; xi1b], -1.0, 1.0)
                l = 1/2*abs(xi1[2]-xi1[1])
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
                    xi = ip.coords[1]
                    xi_s = dot([1/2*(1-xi); 1/2*(1+xi)], xi1)
                    N1 = vec(get_basis(slave_element, xi_s, time))
                    x_s = interpolate(N1, x1) # coordinate in gauss point
                    n_s = interpolate(N1, n1) # normal direction in gauss point
                    t_s = Q'*n_s # tangent direction in gauss point
                    xi_m = project_from_slave_to_master_ad(master_element, x_s, n_s, x2)
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
                if (avg_gap < 1.0e-12) && (std_gap < 1.0e-12)
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
            lan = dot(normals[:,j], la[:,j])
            condition[j] = ForwardDiff.value(lan - gap[1, j])
            is_active[j] = condition[j] > 0
        end

        if problem.properties.iteration == 1 && state == :ACTIVE
            for j in S
                is_active[j] = true
            end
        end

        if problem.properties.iteration == 1 && state == :INACTIVE
            for j in S
                is_active[j] = false
            end
        end

        if problem.properties.print_summary
            @info("Summary of nodes")
            for j in sort(collect(keys(is_active)))
                n = map(ForwardDiff.value, normals[:,j])
                @info("$j, c=$(condition[j]), s=$(is_active[j]), n=$n")
            end
        end

        for j in S

            if is_active[j]
                n = normals[:,j]
                t = Q'*n
                lan = dot(n, la[:,j])
                lat = dot(t, la[:,j])
                C[1,j] -= gap[1, j]
                C[2,j] += lat
            else
                C[:,j] = la[:,j]
            end

        end

        # Apply scaling to contact virtual work and contact constraints
        if problem.properties.use_scaling
            scaling = empty!(problem.properties.scaling_factors)
            for j in S
                isapprox(beta[j], 0.0) && continue
                #is_active[j] || continue
                #haskey(alpha, j) || continue
                #haskey(beta, j) || continue
                #haskey(nadj_nodes, j) || continue
                scaling[j] = alpha[j] / (nadj_nodes[j] * beta[j])
                C[1,j] *= scaling[j]
                fc[:,j] *= scaling[j]
            end
            update!(slave_elements, "contact scaling", time => scaling)
        end

        return vec([fc C])

    end

    # x doesn't mean deformed configuration here
    x = [problem.assembly.u; problem.assembly.la]
    if length(x) == 0
        error("2d autodiff contact problem: initialize problem.assembly.u & la before solution")
    end

    A = ForwardDiff.jacobian(calculate_interface, x)
    b = calculate_interface(x)
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

    f += C1*problem.assembly.la
    g += D*problem.assembly.la

    #=
    if !haskey(problem, "contact force")
        problem.fields["contact force"] = Field(time => f)
    else
        update!(problem.fields["contact force"], time => f)
    end

    fc = problem.fields["contact force"]

    if length(fc) > 1
        # kick in generalized alpha rule for time integration
        alpha = 0.5
        info("Applying Generalized alpha time integration")
        K = (1-alpha)*K
        C1 = (1-alpha)*C1
        f = alpha*fc[end-1].data
    end
    =#

    problem.assembly.K = K
    problem.assembly.C1 = copy(transpose(C1))
    problem.assembly.C2 = C2
    problem.assembly.D = D
    problem.assembly.f = sparse(f)
    problem.assembly.g = sparse(g)

end
