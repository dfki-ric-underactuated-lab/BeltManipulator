"""
    friction_model(q̇, coeff::Tuple)

Accurate friction model considering coulomb and viscous friction 
"""
function friction_model(q̇, coeff::Tuple)
    τ_f = coeff[1]*sign.(q̇) + coeff[2]*q̇
    return τ_f
end   
# # Smooting that function by a sigmoid does not make things better numerically
# function friction_model(q̇, coeff::Tuple; ssc::Real = 100)
#     τ_f = 2coeff[1]/(1 + exp(-ssc*q̇)) -coeff[1] + coeff[2]*q̇
#     return τ_f
# end  

"""
Continuous dynamics of the system using RigidBodyDynamics
"""
function dynamics(mec::RBDMechanism, x::AbstractVector{T1}, u::AbstractVector{T2}) where {T1, T2}
    T = promote_type(T1, T2)
    state = mec.statecache[T]
    res = mec.dyncache[T]

    copyto!(state, x)

    if isnothing(mec.friction)
        dynamics!(res, state, mec.S'*u)
    else
        friction = friction_model.(parent(velocity(state)), mec.friction)
        dynamics!(res, state, mec.S'*u - friction)
    end 
    
    return [res.q̇; res.v̇]
end 
