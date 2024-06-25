"""
    runge_kutta4(mod::Model, x::Vector, u::Vector, h::Float64)
    
4th order Runge-Kutta method with state 'x' and zero-order hold on 'u'.
"""
function runge_kutta4(mod::Model, x::Vector, u::Vector, h::Float64)
    k1 = mod.dynamics(x, u)
    k2 = mod.dynamics(x + 0.5h*k1, u)
    k3 = mod.dynamics(x + 0.5h*k2, u)
    k4 = mod.dynamics(x + h*k3, u)
    return x +  h/6*(k1 + 2k2 + 2k3 + k4)
end

"""
    runge_kutta4(mod::Model, x::Vector, h::Float64)
    
4th order Runge-Kutta method with state 'x' as purely passive dynamics
"""
function runge_kutta4(mod::Model, x::Vector, h::Float64)
    k1 = mod.dynamics(x)
    k2 = mod.dynamics(x + 0.5h*k1)
    k3 = mod.dynamics(x + 0.5h*k2)
    k4 = mod.dynamics(x + h*k3)
    return x +  h/6*(k1 + 2k2 + 2k3 + k4)
end

"""
    linearize_model!(mod::ODEModel, x::Vector{Float64}, u::Vector{Float64}, h::Float64)

Linearizes the dynamics model, discretized by RK4 by automatic differentiation.
"""
function linearize_model!(mod::Model, x::Vector{Float64}, u::Vector{Float64}, h::Float64)
    FD.jacobian!(mod.A, x_ -> runge_kutta4(mod, x_, u, h), x)
    FD.jacobian!(mod.B, u_ -> runge_kutta4(mod, x, u_, h), u)
end

"""
    second_order_derivative!(mod::ODEModel, x::Vector{Float64}, u::Vector{Float64}, h::Float64)

Hessian of the dynamics by autodiff
"""
function second_order_derivative!(mod::ODEModel, x::Vector{Float64}, u::Vector{Float64}, h::Float64)
    FD.jacobian!(mod.Ax, x_ -> vec(FD.jacobian(x__ -> runge_kutta4(mod, x__, u, h), x_)), x)
    FD.jacobian!(mod.Bx, x_ -> vec(FD.jacobian(u_ -> runge_kutta4(mod, x_, u_, h), u)), x)
    # FD.jacobian!(mod.Au, u_ -> vec(FD.jacobian(x_ -> runge_kutta4(mod, x_, u_, h), x)), u)
    FD.jacobian!(mod.Bu, u_ -> vec(FD.jacobian(u__ -> runge_kutta4(mod, x, u__, h), u_)), u)    
end


"""
    simulate!(trj::Trajectory, mod::Model, cp::ControlParameters)

Integrates on 'trj' with its initial value by using a controller 'cp' that potentially depends on the model 'mod'.
"""
function simulate!(trj::Trajectory, mod::Model, cp::ControlParameters)
    t_start = time_ns()
    for k = 1:trj.N-1
        trj.U[k] = return_controls(mod, cp, trj.X[k], k)
        trj.X[k+1] = runge_kutta4(mod, trj.X[k], trj.U[k], trj.times[k+1] - trj.times[k])
    end 
    t_end = time_ns()

    rate = (trj.N - 1)/(t_end - t_start)*1e9
    println("Controller ran at $(rate) Hz")
end


"""
    naive_trajectory!(trj::Trajectory, mec::RBDMechanism, op::iLQRParameters)

Creates a trajectory based on inverse dynamics of the model, by assuming "bang-bang" acceleration in between boundary points.
"""
function naive_trajectory!(trj::Trajectory, mec::RBDMechanism, op::iLQRParameters)

    m = num_positions(mec.rbd_model)

    # interpolation of joint coordinates 
    Δθ = op.xf[1:m] - op.x0[1:m]

    # assuming "bang-bang" acceleration vector
    a = 4Δθ/trj.times[end]^2

    # fake_state is used to get the acceleration vector as SegmentedVector, which is not sufficiently documented to be created alone
    fake_state = MechanismState(mec.rbd_model, zeros(m), a)
    fake_state_ = MechanismState(mec.rbd_model, zeros(m), -a)
    
    for k = 1:trj.N-1
        t = trj.times[k]
        if t <= trj.times[end]/2
            pos = 0.5*t^2*a + op.x0[1:m]
            vel = t*a
            current_state = MechanismState(mec.rbd_model, pos, vel)
            joint_torques = inverse_dynamics(current_state, fake_state.v)
            if !isnothing(mec.friction)
                joint_torques += friction_model.(vel, mec.friction)
            end
            trj.U[k] = mec.S_'*joint_torques
            trj.X[k] = [pos; vel]
        else
            pos = op.x0[1:m] - a/4*trj.times[end]^2 + a*(t*trj.times[end] - 0.5t^2)
            vel = a*(trj.times[end]-t)
            current_state = MechanismState(mec.rbd_model, pos, vel) 
            joint_torques = inverse_dynamics(current_state, fake_state_.v)
            if !isnothing(mec.friction)
                joint_torques += friction_model.(vel, mec.friction)
            end
            trj.U[k] = mec.S_'*joint_torques
            trj.X[k] = [pos; vel]
        end
    end
    trj.X[end] = op.xf
end

"""
    torques_from_trajectory!(trj::Trajectory, man::Manipulator)

Computes from a known joint-space trajectory the torques in motor space
"""
function torques_from_trajectory!(trj::Trajectory, mec::RBDMechanism)

    m = num_positions(mec.rbd_model)
    fake_state = MechanismState(mec.rbd_model, zeros(m), zeros(m))
    accs = numerical_time_derivative([trj.X[k][m+1:end] for k = 1:trj.N], trj.times)

    for k = 1:trj.N - 1
        current_state = MechanismState(mec.rbd_model, trj.X[k][1:m], trj.X[k][m+1:end])
        set_velocity!(fake_state, accs[k])
        joint_torques = inverse_dynamics(current_state, fake_state.v)
        if !isnothing(mec.friction)
            joint_torques += friction_model.(trj.X[k][m+1:end], mec.friction)
        end
        trj.U[k] = mec.S_'*joint_torques
    end
end