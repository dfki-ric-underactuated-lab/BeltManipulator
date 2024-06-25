
"""
    trajectory_cost(trj::Trajectory, op::iLQRParameters)

Returns the additive cost and final cost at the goal state
"""
function trajectory_cost(trj::Trajectory, op::iLQRParameters)
    J = 0.
    # stage cost
    for k = 1:(trj.N-1)
        J += 0.5*((trj.X[k] - op.xf)'*op.Q*(trj.X[k] - op.xf)) + 0.5*(trj.U[k]'*op.R*trj.U[k])
    end
    # terminal cost
    J += 0.5*((trj.X[end] - op.xf)'*op.Qf*(trj.X[end] - op.xf))
    return J
end

"""
    backward_pass!(trj::Trajectory, mod::Model, op::iLQRParameters, fpf::Bool = false)

Riccati sweep for a given trajectory, based on the weight matrices. The optimization parameters 'op' are altered with updated gains and feedforward values.
"""
function backward_pass!(op::iLQRParameters, trj::Trajectory, model::Model)

    # storing two terms for expected cost change
    ΔVu = 0.
    ΔVuu = 0.

    # optimal terminal cost to go
    p = op.Qf*(trj.X[end] - op.xf)
    P = op.Qf

    k = trj.N-1
    reg_limit = 0
    while k > 0
        # getting A and B matrices
        linearize_model!(model, trj.X[k], trj.U[k], trj.times[k+1] - trj.times[k])
        
        # computing linearized action-value function
        Qx = op.Q*(trj.X[k] - op.xf) + model.A'*p
        Qu = op.R*trj.U[k] + model.B'*p
        
        if op.second_order
            second_order_derivative!(model, trj.X[k], trj.U[k], trj.times[k+1] - trj.times[k])
            Qxx = op.Q + model.A'*P*model.A + kron(p', I(model.m))*commutation_matrix(model.m, model.m)*model.Ax
            Qux = model.B'*P*model.A + kron(p', I(model.n))*commutation_matrix(model.m, model.n)*model.Bx
            Quu = op.R + model.B'*P*model.B + kron(p', I(model.n))*commutation_matrix(model.m, model.n)*model.Bu
        else
            Qxx = op.Q + model.A'*P*model.A
            Qux = model.B'*P*model.A
            Quu = op.R + model.B'*P*model.B 
        end

        # regularization step
        if !iszero(op.β)
            if op.regularization == :state
                IB = model.B'*op.β*I*model.B
                IA = model.B'*op.β*I*model.A
                Qux = Qux .+ IA
            elseif op.regularization == :control
                IB = op.β*I(model.n)
            end
            Quu = Quu .+ IB
        end

        # repeating backward pass, if necessary
        # Quu_fact = chol!(Quu)::Cholesky
        # if !isposdef(Quu_fact)
        if !isinvertible(Quu) && reg_limit < op.reg_max
            @warn "Backward pass repeated at time step $k in iteration $(op.count)"
            # moderate regularization increase
            op.β += op.dβ
            k = trj.N-1
            reg_limit += 1
            ΔVu = 0.; ΔVuu = 0.
            continue
        elseif reg_limit == op.reg_max
            @warn "Backward pass failed in iteration $(op.count)"
            reg_limit = 0
            op.β = 0.
            break
        end
        
        op.d[k] = -Quu\Qu
        op.K[k] = -Quu\Qux

        # expected change in cost-to-go
        ΔVu += op.d[k]'*Qu
        ΔVuu += 0.5*op.d[k]'*Quu*op.d[k]

        # updating p and P
        p = Qx + op.K[k]'*Quu*op.d[k] + op.K[k]'*Qu + Qux'*op.d[k]
        P = Qxx + op.K[k]'*Quu*op.K[k] + op.K[k]'*Qux + Qux'*op.K[k]

        k += -1
    end

    return ΔVu, ΔVuu
end

"""
    forward_pass!(trj::Trajectory, mod::Model, J::Float64, ΔVu::Float64, ΔVuu::Float64, op::iLQRParameters)

Applies the feedback and feedforward values onto the controls and makes forward rollout
"""
function forward_pass!(trj::Trajectory, model::Model, J::Float64, ΔVu::Float64, ΔVuu::Float64, op::iLQRParameters)

    # copying current state before starting line search
    copyto!(op.trj_new.X[1], trj.X[1])

    # forward rollout and line search
    α = 1.
    line_search_count = 0
    Jn = 0.
    while true
        if op.fail_max == op.fails
            @warn "Maximum fails of $(op.fail_max) reached - stopping now"
            break
        end

        try
            for k = 1:trj.N-1
                op.trj_new.U[k] = trj.U[k] + op.K[k]*(op.trj_new.X[k] - trj.X[k]) + α*op.d[k]
                op.trj_new.X[k+1] = runge_kutta4(model, op.trj_new.X[k], op.trj_new.U[k], trj.times[k+1] - trj.times[k])
            end
        catch
            @warn "Rollout failed at iteration $(op.count)!"
            # changing line-search parameter
            α = op.a_red*α
            op.fails += 1
            continue
        end

        Jn = trajectory_cost(op.trj_new, op)
 
        if op.z_low < -(J - Jn)/(α*ΔVu + α^2*ΔVuu) < op.z_upp
            copy_traj!(trj, op.trj_new)
            op.β = 0
            op.fails = 0
            break
        elseif line_search_count == op.ls_max
            @warn "Line-search failed at iteration $(op.count)!"
            # moderate regularization increase to repeat backward pass
            op.β += op.dβ
            op.count += -1
            op.fails += 1
            break
        else
            α = op.a_red*α
            line_search_count += 1
        end
                 
    end
    return Jn
end 

"""
    iLQR!(trj::Trajectory, model::Model, op::iLQRParameters)
    
Iteratively changes the trajectory 'trj' by means of the iLQR. SimulationParameters have to be supplied for using the generalized-α method in the forward pass.
""" 
function iLQR!(trj::Trajectory, model::Model, op::iLQRParameters)   
    # performing initial rollout with given controls
    copyto!(trj.X[1], op.x0)

    for k = 1:trj.N-1
        trj.X[k+1] = runge_kutta4(model, trj.X[k], trj.U[k], trj.times[k+1] - trj.times[k])
    end
    
    ΔVu, ΔVuu = backward_pass!(op, trj, model)
    J = trajectory_cost(trj, op)
    op.count = 1

    while maximum(abs.(collect(Iterators.flatten(op.d)))) > op.tol
        op.count += 1
        J = forward_pass!(trj, model, J, ΔVu, ΔVuu, op)

        ΔVu, ΔVuu = backward_pass!(op, trj, model)

        if op.count >= op.i_max || op.fail_max == op.fails
            break
        end
    end
    
    println("Took $(op.count) iterations of max $(op.i_max)")

end