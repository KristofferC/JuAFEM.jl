using JuAFEM, SparseArrays, BlockArrays

import DifferentialEquations

grid = generate_grid(Quadrilateral, (60, 60), Vec{2}((0.0,0.0)), Vec{2}((2.5,2.5)))
addnodeset!(grid, "ground", x-> x[2] == -0 && x[1] == -0)
dim = 2
Δt = 0.1
T = 1000
ip = Lagrange{dim, RefCube, 1}()
qr = QuadratureRule{dim, RefCube}(2)
cellvalues = CellScalarValues(qr, ip);

dh = DofHandler(grid)
push!(dh, :ϕₘ, 1)
push!(dh, :ϕₑ, 1)
push!(dh, :s, 1)
close!(dh);

K = create_sparsity_pattern(dh)
M = create_sparsity_pattern(dh);

Base.@kwdef struct FHNParameters
    a::Float64 = 0.1
    b::Float64 = 0.5
    c::Float64 = 1.0
    d::Float64 = 0.0
    e::Float64 = 0.01
end;

function κₑ(x)
    return SymmetricTensor{2,2,Float64}((3.5e-5, 0, 2.5e-5))
end;

function κᵢ(x)
    return SymmetricTensor{2,2,Float64}((4.5e-5, 0, 2.0e-6))
end;

function Cₘ(x)
    return 1.0
end;

function χ(x)
    return 1.0
end;

ch = ConstraintHandler(dh)
∂Ω = getnodeset(grid, "ground")
dbc = Dirichlet(:ϕₑ, ∂Ω, (x, t) -> 0)
add!(ch, dbc)
close!(ch)
update!(ch, 0.0);

function doassemble_linear!(cellvalues::CellScalarValues{dim}, K::SparseMatrixCSC, M::SparseMatrixCSC, dh::DofHandler; params::FHNParameters = FHNParameters()) where {dim}
    n_ϕₘ = getnbasefunctions(cellvalues)
    n_ϕₑ = getnbasefunctions(cellvalues)
    n_s = getnbasefunctions(cellvalues)
    ntotal = n_ϕₘ + n_ϕₑ + n_s
    n_basefuncs = getnbasefunctions(cellvalues)
    #We use PseudoBlockArrays to write into the right places of Ke
    Ke = PseudoBlockArray(zeros(ntotal, ntotal), [n_ϕₘ, n_ϕₑ, n_s], [n_ϕₘ, n_ϕₑ, n_s])
    Me = PseudoBlockArray(zeros(ntotal, ntotal), [n_ϕₘ, n_ϕₑ, n_s], [n_ϕₘ, n_ϕₑ, n_s])

    assembler_K = start_assemble(K)
    assembler_M = start_assemble(M)

    #Here the block indices of the variables are defined.
    ϕₘ▄, ϕₑ▄, s▄ = 1, 2, 3

    #Now we iterate over all cells of the grid
    @inbounds for cell in CellIterator(dh)
        fill!(Ke, 0)
        fill!(Me, 0)
        #get the coordinates of the current cell
        coords = getcoordinates(cell)

        JuAFEM.reinit!(cellvalues, cell)
        #loop over all Gauss points
        for q_point in 1:getnquadpoints(cellvalues)
            #get the spatial coordinates of the current gauss point
            coords_qp = spatial_coordinate(cellvalues, q_point, coords)
            #based on the gauss point coordinates, we get the spatial dependent
            #material parameters
            κₑ_loc = κₑ(coords_qp)
            κᵢ_loc = κᵢ(coords_qp)
            Cₘ_loc = Cₘ(coords_qp)
            χ_loc = χ(coords_qp)
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:n_basefuncs
                Nᵢ = shape_value(cellvalues, q_point, i)
                ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
                for j in 1:n_basefuncs
                    Nⱼ = shape_value(cellvalues, q_point, j)
                    ∇Nⱼ = shape_gradient(cellvalues, q_point, j)
                    #diffusion parts
                    Ke[BlockIndex((ϕₘ▄,ϕₘ▄),(i,j))] -= ((κᵢ_loc ⋅ ∇Nᵢ) ⋅ ∇Nⱼ) * dΩ
                    Ke[BlockIndex((ϕₘ▄,ϕₑ▄),(i,j))] -= ((κᵢ_loc ⋅ ∇Nᵢ) ⋅ ∇Nⱼ) * dΩ
                    Ke[BlockIndex((ϕₑ▄,ϕₘ▄),(i,j))] -= ((κᵢ_loc ⋅ ∇Nᵢ) ⋅ ∇Nⱼ) * dΩ
                    Ke[BlockIndex((ϕₑ▄,ϕₑ▄),(i,j))] -= (((κₑ_loc + κᵢ_loc) ⋅ ∇Nᵢ) ⋅ ∇Nⱼ) * dΩ
                    #linear reaction parts
                    Ke[BlockIndex((ϕₘ▄,ϕₘ▄),(i,j))] -= params.a * Nᵢ * Nⱼ * dΩ
                    Ke[BlockIndex((ϕₘ▄,s▄),(i,j))]  -= Nᵢ * Nⱼ * dΩ
                    Ke[BlockIndex((s▄,ϕₘ▄),(i,j))]  += params.e * params.b * Nᵢ * Nⱼ * dΩ
                    Ke[BlockIndex((s▄,s▄),(i,j))]   -= params.e * params.c * Nᵢ * Nⱼ * dΩ
                    #mass matrices
                    Me[BlockIndex((ϕₘ▄,ϕₘ▄),(i,j))] += Cₘ_loc * χ_loc * Nᵢ * Nⱼ * dΩ
                    Me[BlockIndex((s▄,s▄),(i,j))]   += Nᵢ * Nⱼ * dΩ
                end
            end
        end

        assemble!(assembler_K, celldofs(cell), Ke)
        assemble!(assembler_M, celldofs(cell), Me)
    end
    return K, M
end;

function apply_nonlinear!(du, u, p, t)
    dh = p[2]
    ch = p[3]
    params = p[4]
    ip = p[5]
    qr = p[6]
    cellvalues = p[7]
    n_basefuncs = getnquadpoints(cellvalues)

    for cell in CellIterator(dh)
        JuAFEM.reinit!(cellvalues, cell)
        _celldofs = celldofs(cell)
        ϕₘ_celldofs = _celldofs[dof_range(dh, :ϕₘ)]
        s_celldofs = _celldofs[dof_range(dh, :s)]
        ϕₘe = u[ϕₘ_celldofs]
        se = u[s_celldofs]
        coords = getcoordinates(cell)
        for q_point in 1:getnquadpoints(cellvalues)
            x_qp = spatial_coordinate(cellvalues, q_point, coords)
            χ_loc = χ(x_qp)
            dΩ = getdetJdV(cellvalues, q_point)
            val = function_value(cellvalues, q_point, ϕₘe)
            nl_contrib = - val^3 + (1+params.a)*val^2
            for j in 1:n_basefuncs
                Nⱼ = shape_value(cellvalues, q_point, j)
                du[ϕₘ_celldofs[j]] += χ_loc * nl_contrib * Nⱼ * dΩ
                du[s_celldofs[j]]  -= params.e * params.d * Nⱼ * dΩ
            end
        end
    end
    apply_zero!(du, ch)
end;

K, M = doassemble_linear!(cellvalues, K, M, dh);

apply!(K, ch)
apply!(M, ch);

function bidomain!(du,u,p,t)
    du .= K * u
    apply_nonlinear!(du, u, p, t)
end;

u₀ = zeros(ndofs(dh));
for cell in CellIterator(dh)
    _celldofs = celldofs(cell)
    ϕₘ_celldofs = _celldofs[dof_range(dh, :ϕₘ)]
    s_celldofs = _celldofs[dof_range(dh, :s)]
    for (i, coordinate) in enumerate(getcoordinates(cell))
        if coordinate[2] >= 1.25
            u₀[s_celldofs[i]] = 0.1
        end
        if coordinate[1] <= 1.25 && coordinate[2] <= 1.25
            u₀[ϕₘ_celldofs[i]] = 1.0
        end
    end
end

jac_sparsity = sparse(K)

f = DifferentialEquations.ODEFunction(bidomain!,mass_matrix=M;jac_prototype=jac_sparsity)
p = [K, dh, ch, FHNParameters(), ip, qr, cellvalues]
prob_mm = DifferentialEquations.ODEProblem(f,u₀,(0.0,T),p)

sol = DifferentialEquations.solve(prob_mm,DifferentialEquations.QNDF(),reltol=1e-3,abstol=1e-4, progress=true, progress_steps = 1, adaptive=true, dt=Δt);

pvd = paraview_collection("bidomain.pvd");

for (solution,t) in zip(sol.u, sol.t)
    #compress=false flag because otherwise each vtk file will be stored in memory
    vtk_grid("bidomain-$t.vtu", dh; compress=false) do vtk
        vtk_point_data(vtk,dh,solution)
        vtk_save(vtk)
        pvd[t] = vtk
    end
end

vtk_save(pvd);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

