"""
Abstraction and iterable for dynamic models and their linearizations
"""
abstract type Model end;

mutable struct ODEModel <: Model
    const m::Int                                        # state-space variables                 
    const n::Int                                        # controls

    dynamics::Function
    
    # state variables
    x::Vector{Float64}
    u::Vector{Float64}
    
    # Linearized dynamics in state-space form
    A::Matrix{Float64}
    B::Matrix{Float64}

    # second order derivatives in vector form
    Ax::Matrix{Float64}
    Bx::Matrix{Float64}
    # Au::Matrix{Float64}
    Bu::Matrix{Float64}

    """
        ODEModel(m::Int, n::Int, dynamics::Function)

    Constructor for ODE models that are solely defined by the state-space dynamics 'dynamics'. 'm' is number of states and 'n' of controls.
    """
    function ODEModel(m::Int, n::Int, dynamics::Function)

        A = fill(NaN, (m, m))
        B = fill(NaN, (m, n))
        Ax = fill(NaN, (m*m, m))
        Bx = fill(NaN, (m*n, m))
        # Au = fill(NaN, (m*m, n))
        Bu = fill(NaN, (m*n, n))


        return new(m, n, dynamics, fill(NaN, m), fill(NaN, n), A, B, Ax, Bx, Bu)
    end
end

function ODEModel(n::Int, mec::RBDMechanism)
    dyn = (x, u) -> dynamics(mec, x, u)
    ODEModel(2num_positions(mec.rbd_model), n, dyn)
end