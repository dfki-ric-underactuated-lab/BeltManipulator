"""
Returns the interval values for states and controls at time-step 'k'
"""
function return_indices(i::Int, m::Int, n::Int)
    ri = (i - 1)*(m + n)
    xs = ri+n+1:ri+m+n
    us = ri+1:ri+n
    return (xs, us)
end

"""
Updates the QP at the time-step k with the pre-computed linearized terms
"""
function setup_solver!(cp::MPCParameters, mod::Model, x::Vector{Float64}, k::Int, r::Int)
    n, m = cp.n, cp.m

    # getting reference points from controller
    xref = cp.Tref.X[r]
    uref = cp.Tref.U[r]

    # linearizing at the current point and taking previous controls, if available
    if k == 1
        u = zeros(n)
    else
        u = cp.osqp_sol[1:n]
    end
    linearize_model!(mod, x, u, cp.Tref.times[k+1] - cp.Tref.times[k])
    k_vec = runge_kutta4(mod, x, u, cp.Tref.times[k+1] - cp.Tref.times[k])
    k_vec += -mod.A*x + -mod.B*u


    # changing first rows and columns
    cp.A[1:m,1:n] .= mod.B

    for i = 1:cp.hx-1
        xs, us = return_indices(i, m, n)
        # linear constraints - linearized dynamics at current position and constant vector
        cp.A[i*m+1:(i+1)*m, (i-1)*(m+n)+n+1:i*(m+n)] .= mod.A
        cp.A[i*m+1:(i+1)*m, i*(m+n)+1:i*(m+n)+n] .= mod.B

        # linear cost terms are updated
        cp.q[us] .= -cp.R*uref
        cp.q[xs] .= -cp.Q*xref
    end

    # updating first vector entry of linearized dynamics
    cp.lb[1:m] .= -mod.A*x - k_vec
    cp.ub[1:m] .= -mod.A*x - k_vec
    # updating remaining part (H-1)*n in bounds vectors with constant values from linearized dynamics
    cp.lb[m+1:cp.hx*m] .= repeat(-k_vec, cp.hx-1)
    cp.ub[m+1:cp.hx*m] .= repeat(-k_vec, cp.hx-1)

    # overwriting last entry of q-vector
    cp.q[end-m+1:end] .= -cp.Qf*xref

    # writing into OSQP
    OSQP.update!(cp.solver, q=cp.q, Ax=sparse(cp.A).nzval, l=cp.lb, u=cp.ub)
end

"""
    return_controls(mod::Model, cp::MPC1Parameters, x::Vector{Float64}, k::Int)

Returns the control values by solving the QP at time-step 'k' and state 'x'
"""
function return_controls(mod::Model, cp::MPCParameters, x::Vector{Float64}, k::Int)
    
    if cp.count == cp.hu
        # finding correct reference point on trajectory
        r = min(k + cp.hx, cp.Tref.N - 1)
        # clipping k that we don't try to control from the 'future'
        k = min(k, cp.Tref.N - 1)

        setup_solver!(cp, mod, x, k, r)
        res = OSQP.solve!(cp.solver)
        cp.osqp_sol = res.x
        cp.count = 1
        return cp.osqp_sol[1:cp.n]
    else
        cp.count += 1
        return cp.osqp_sol[(cp.m+cp.n)*(cp.count-1)+1:(cp.m+cp.n)*(cp.count-1)+cp.n]
    end
end

"""
    call_controller(mod::Model, cp::ControlParameters, x::Vector{Float64})

Returns motor torques from the current state 'x' for a given controller, defined under ControlParameters that contains a reference Trajectory
"""
function call_controller(mod::Model, cp::ControlParameters, x::Vector{Float64})
    # storing start time of the controller
    if cp.t0 === nothing
        cp.t0 = Int(time_ns())
    end

    # getting the correct reference
    return return_controls(mod, cp, x, searchsortedlast(cp.Tref.times, (Int(time_ns()) - cp.t0)/1e9))
end

"""
    return_controls(mod::Model, cp::FeedForward, x::Vector, k::Int)

Ignores model and state and simply returns the reference controls
"""
function return_controls(mod::Model, cp::FeedForward, x::Vector, k::Int)
    k = min(k, cp.Tref.N - 1)
    return cp.Tref.U[k]
end

"""
    return_controls(mod::Model, cp::NoControl, x::Vector, k::Int)

No control at all
"""
function return_controls(mod::Model, cp::NoControl, x::Vector, k::Int)
    return zeros(mod.n)
end
