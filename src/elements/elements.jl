# Classic structural elements
include("CALFEM_API.jl")
include("CALFEM_docs.jl")
include("spring.jl")
include("bar.jl")

include("solid_elements.jl")
include("heat_elements.jl")

create_initial_vars(inits) = Expr(:block, [:($sym = zeros($size)) for (sym, size) in inits]...)

"""
This function generates the body expression for a function that
computes the so called stiffness matrix and force vector
for a finite element.

The finite element should be an instance of `FElement`.
"""
function gen_ke_fe_body(ele)
    quote
        ndofs = $(ele.nnodes) * $(ele.dofs_per_node)
        nnodes = $(ele.nnodes)
        ndim = $(get_ndim(ele))

        # TODO; need to fix the eq size for different problems
        #length(eq) == ndim || throw(ArgumentError("length of eq must be $ndim"))
        int_order > 0 || throw(ArgumentError("integration order must be > 0"))

        # Default buffers
        Ke = zeros(ndofs, ndofs)
        fe = zeros(ndofs)
        GRAD_KERNEL = zeros(ndofs, ndofs)
        SOURCE_KERNEL = zeros(ndofs)
        dNdx = zeros(ndim, nnodes)
        dNdξ = zeros(ndim, nnodes)
        J = zeros(ndim, ndim)
        N = zeros(nnodes)

        # Create the element requested buffers
        $(create_initial_vars(ele.inits))

        # Check if we need to compute the RHS
        compute_RHS = (eq != zeros(ndim))

        # Create the gauss rule of the requested order on the
        # elements reference shape
        qr = get_gaussrule($(ele.shape), int_order)

        for (ξ, w) in zip(qr.points, qr.weights)
            evaluate_N!($(ele.shape_funcs), N, ξ)
            evaluate_dN!($(ele.shape_funcs), dNdξ, ξ)
            @into! J = dNdξ * x
            Jinv = inv_spec(J)
            @into! dNdx = Jinv * dNdξ
            dV = det_spec(J) * w

            ##############################
            # Call the elements LHS kernel
            $(ele.grad_kernel())
            ##############################

            @devec Ke[:, :] += GRAD_KERNEL .* dV

            if compute_RHS
                ##############################
                # Call the elements RHS kernel
                $(ele.source_kernel())
                ##############################

                @devec fe += SOURCE_KERNEL .* dV
            end
        end
        return Ke, fe
    end
end

for fem in [S_S_1, S_S_2, S_T_1, S_C_1, # Solid elements
            H_S_1, H_S_2, H_T_1, H_C_1]
    f_name = Symbol(string(fem.name) * "e")
    if get_ndim(fem) == 2
        @eval function $(f_name)(x::Matrix, D::Matrix, t::Number,
                                  eq::VecOrMat=zeros($(get_ndim(fem))),
                                  int_order::Int=$(fem.default_intorder))
            $(gen_ke_fe_body(fem))
        end
    elseif get_ndim(fem) == 3
        @eval function $(f_name)(x::Matrix, D::Matrix,
                                  eq::VecOrMat=zeros($(get_ndim(fem))),
                                  int_order::Int=$(fem.default_intorder))
            $(gen_ke_fe_body(fem))
        end
    end
end


"""
This function generates the body expression for a function that
computes the flux (stress) and its conjugate(strain)

The finite element should be an instance of `FElement`.
"""
function gen_s_body(ele)
    quote
        ndofs = $(ele.nnodes) * $(ele.dofs_per_node)
        nnodes = $(ele.nnodes)
        ndim = $(get_ndim(ele))
        fluxlen = $(get_nflux(ele))

        #length(eq) == ndim || throw(ArgumentError("length of eq must be $ndim"))
        int_order > 0 || throw(ArgumentError("integration order must be > 0"))

        # Default buffers
        FLUX_KERNEL = zeros(fluxlen)
        CONJ_KERNEL = zeros(fluxlen)
        dNdx = zeros(ndim, nnodes)
        dNdξ = zeros(ndim, nnodes)
        J = zeros(ndim, ndim)
        N = zeros(nnodes)

        # Create the element requested buffers
        $(create_initial_vars(ele.inits))

        # Create the gauss rule of the requested order on the
        # elements reference shape
        qr = get_gaussrule($(ele.shape), int_order)

        n_gps = length(qr.points)
        FLUXES = zeros(n_gps, fluxlen)
        CONJS = zeros(n_gps, fluxlen)
        points = zeros(n_gps, ndim)

        for (i, ξ) in enumerate(qr.points)
            evaluate_N!($(ele.shape_funcs), N, ξ)
            evaluate_dN!($(ele.shape_funcs), dNdξ, ξ)
            @into! J = dNdξ * x
            Jinv = inv_spec(J)
            @into! dNdx = Jinv * dNdξ

            ##############################
            # Call the elements flux kernel
            $(ele.flux_kernel())
            ##############################

            FLUXES[i, :] = FLUX_KERNEL
            CONJS[i, :] = CONJ_KERNEL
            points[i, :] = x' * N
        end
        return FLUXES, CONJS, points
    end
end

for fem in [S_S_1, S_S_2, S_T_1, S_C_1, # Solid elements
            H_S_1, H_S_2, H_T_1, H_C_1]
    f_name = Symbol(string(fem.name) * "s")
    if get_ndim(fem) == 2
        @eval function $(f_name)(x::Matrix, D::Matrix, t::Number,
                                 ed::VecOrMat,
                                 int_order::Int=$(fem.default_intorder))
            $(gen_s_body(fem))
        end
    elseif get_ndim(fem) == 3
        @eval function $(f_name)(x::Matrix, D::Matrix,
                                 ed::VecOrMat,
                                 int_order::Int=$(fem.default_intorder))
            $(gen_s_body(fem))
        end
    end
end
