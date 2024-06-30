mutable struct RBDMechanism 

    rbd_model::RigidBodyDynamics.Mechanism     
    statecache::RigidBodyDynamics.StateCache                               
    dyncache::RigidBodyDynamics.DynamicsResultCache    
    
    # visualization
    vis::Union{Nothing, Visualizer};
    mvis::Union{Nothing, MechanismVisualizer};
    
    friction::Union{Nothing, Vector{Tuple{Float64, Float64}}} 

    # structure matrix, if mechanism is belt-driven
    S::Matrix{Float64}
    S_::Matrix{Float64}

    function RBDMechanism(urdf::String; S = nothing, friction::Union{Nothing, Vector{Tuple{Float64, Float64}}}  = nothing, visualizer::Bool = true)

        rbd_model = parse_urdf(urdf)
        statecache = StateCache(rbd_model)
        dyncache = DynamicsResultCache(rbd_model)        

        # creating visualizer object and rendering it
        if visualizer
            vis = Visualizer()
            render(vis)
            # create mechanism visualizer object
            mvis = MechanismVisualizer(rbd_model, URDFVisuals(urdf), vis) 
        else
            vis = nothing; mvis = nothing
        end

        # checking, if structure matrix and friction vector is of correct dimension
        dof = num_positions(rbd_model)
        if isnothing(S)
            S = diagm(ones(dof))
        else
            @assert size(S)[2] == dof "Mismatch in RBD model with $(dof) DOF and S-matrix of size $(size(S))"
        end

        if !isnothing(friction)
            @assert length(friction) == num_positions(rbd_model) "Friction vector has to contain $(num_positions(rbd_model)) tuples for (static, viscous)"
        end
        
        return new(rbd_model, statecache, dyncache, vis, mvis, friction, S, pinv(S))
    end
end
