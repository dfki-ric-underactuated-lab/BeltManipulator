"""
All controllers are accomodated under this type
"""
abstract type ControlParameters end

"""
Contains the MPC parameters along with the reference trajectories
"""
mutable struct MPCParameters <: ControlParameters
    const Tref::Trajectory                  # reference trajectory

    const m::Int                            # number of states
    const n::Int                            # number of controls

    Q::AbstractArray                        # cost matrix for states
    R::AbstractArray                        # cost matrix for controls
    Qf::AbstractArray                       # cost matrix for final state

    t0::Union{Int, Nothing}                 # time stamp when controller is called first time

    hx::Int                                 # horizon length of the controller
    hu::Int                                 # interval in which the controller is active
    count::Int                              # counter when to recompute the solution

    solver                                  # solver for the QP problem
    osqp_sol::Vector{Float64}               # solution of the solver
    tol::Float64                            # tolerance for OSQP
    verbose::Bool                           # verbose mode for OSQP

    # intermediate arrays that are fed to the solver
    P::SparseMatrixCSC{Float64,Int}
    q::Vector{Float64}
    A::Matrix{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}

    function MPCParameters(mod::Model; hx::Int, Q::AbstractArray, R::AbstractArray, Qf::AbstractArray, Tref::Trajectory, tol::Float64, hu::Int = 1, xl::Vector{Float64} = -Inf*ones(mod.m), xu::Vector{Float64} = Inf*ones(mod.m), τ_limit::Real = Inf, verbose::Bool = true)
        
        n, m = mod.n, mod.m
        @assert length(xl) == m && length(xu) == m "Upper or lower state bounds are incorrect"
        if size(Qf) != size(Q)
            throw(error("Final and intermediate cost matrix dimensions do not match. Qf ∈ ℝ^{$(size(Qf)[1]) × $(size(Qf)[2])} but Q ∈ ℝ^{$(size(Q)[1]) × $(size(Q)[2])}"))
        end 
        
        # determining overall size of linear constraint matrix
        np = hx*(n + m)
        if τ_limit == Inf
            nd = hx*n
        else
            nd = np
        end
        if !all(xl .== -Inf) && !all(xu .== Inf)
            nd += hx*n
        end

        q = zeros(np)
        # building pattern of initial constraint matrix
        A = zeros(nd,np)
        A[1:m,1:n] .= fill(1, (m,n))        # first B matrix
        A[1:m,n+1:m+n] .= -I(m)             # first identity matrix
        for i = 1:hx-1
            A[i*m+1:(i+1)*m, (i-1)*(m+n)+n+1:i*(m+n)] .= fill(1, (m,m))
            A[i*m+1:(i+1)*m, i*(m+n)+1:i*(m+n)+n] .= fill(1, (m,n))
            A[i*m+1:(i+1)*m, i*(m+n)+n+1:(i+1)*(m+n)] .= -I(m)
        end

        lb = zeros(nd)
        ub = zeros(nd)

        # adding torque and state constraints in linear matrix A and in lb and ub vectors
        if τ_limit && !all(xl .== -Inf) && !all(xu .== Inf)
            offset = hx*m
            A[offset+1:end,:] .= I(hx*(n+m))
            lb[offset+1:end] .= repeat([-man.τ_limit*ones(n); xl], hx)
            ub[offset+1:end] .= repeat([man.τ_limit*ones(n); xu], hx)
        elseif τ_limit
            offset = hx*m
            [A[offset+n*(i-1)+1:offset+n*i, (m+n)*(i-1)+1:(m+n)*(i-1)+n] .= I(n) for i = 1:hx]
            lb[offset+1:end] .= repeat(-man.τ_limit*ones(n), hx)
            ub[offset+1:end] .= repeat(man.τ_limit*ones(n), hx)
        end

        # initializing cost matrix for horizon length, since we don't change the weights
        P = spzeros(np,np)
        for i = 1:hx-1
            xs, us = return_indices(i, m, n)
            P[xs, xs] .= Q
            P[us, us] .= R
        end

        # changing last rows and columns of cost matrix
        P[end-(n+m)+1:end-m, end-(n+m)+1:end-m] .= R
        P[end-m+1:end, end-m+1:end] .= Qf

        # creating model and writing first parameters to it
        solver = OSQP.Model()
        OSQP.setup!(solver, P=P, q=q, A=SparseMatrixCSC(A), l=lb, u=ub, verbose=verbose, eps_rel=tol, polish=1)

        return new(deepcopy(Tref), m, n, Q, R, Qf, nothing, hx, hu, hu, solver, fill(NaN, np), tol, verbose, P, q, A, lb, ub)
    end
end

"""
Controller based on time-varying linear quadratic regulator that is based on precalculated gains
"""
mutable struct TVLQRParameters <: ControlParameters
    const Tref::Trajectory                  # reference trajectory

    const m::Int                            # number of states
    const n::Int                            # number of controls

    t0::Union{Int, Nothing}                 # time stamp when controller is called first time

    A::Vector{Matrix{Float64}}              # vector of state linearization
    B::Vector{Matrix{Float64}}              # vector of control linearization

    P::Vector{Matrix{Float64}}              # vector of cost-to-go interpolation
    K::Vector{Matrix{Float64}}              # vector of feedback matrices
    
    function TVLQRParameters(mod::Model; Q::AbstractArray, R::AbstractArray, Qf::AbstractArray, Tref::Trajectory)
        
        # linearizing the entire trajectory
        A = Vector{Matrix{Float64}}(fill(fill(NaN, size(mod.A)), Tref.N-1))
        B = Vector{Matrix{Float64}}(fill(fill(NaN, size(mod.B)), Tref.N-1))
        # kref = Vector{Vector{Float64}}(fill(fill(NaN, man.n), trj.N-1))

        for k = 1:Tref.N-1
            linearize_model!(man, Tref.X[k], Tref.U[k], Tref.times[k+1] - Tref.times[k])
            A[k] = deepcopy(mod.A)
            B[k] = deepcopy(mod.B)
            # kref[k] = rk4((x_, u_) -> dynamics(man, x_, u_), trj.X[k], trj.U[k], trj.times[k+1] - trj.times[k])
            # kref[k] += -man.A*trj.X[k] + -man.B*trj.U[k]
        end

        # computing the gains of the controller
        K = [fill(NaN, (mod.m, mod.n)) for _ = 1:Tref.N-1]
        P = [fill(NaN, (mod.m, mod.m)) for _ = 1:Tref.N]

        P[end] = Qf
        for k = reverse(1:Tref.N-1)
            K[k] = (R + B[k]'P[k+1]*B[k])\(B[k]'P[k+1]*A[k])
            P[k] = Q + A[k]'P[k+1]*A[k] - A[k]'P[k+1]*B[k]*K[k]
        end

        return new(deepcopy(Tref), mod.m, mod.n, nothing, A, B, P, K)
    end
    
end

"""
No controller, just feeds reference controls to the system
"""
mutable struct FeedForward <: ControlParameters
    const Tref::Trajectory                  # reference trajectory

    t0::Union{Int, Nothing}                 # time stamp when controller is called first time

    function FeedForward(Tref::Trajectory)
        return new(deepcopy(Tref), nothing)
    end
end

"""
No controller, used for passive dynamics
"""
mutable struct NoControl <: ControlParameters
    const Tref::Trajectory                  # reference trajectory
    
    t0::Union{Int, Nothing}                 # time stamp when controller is called first time

    function NoControl()
        return new(Trajectory(0,0,0,0:0), nothing)
    end
end