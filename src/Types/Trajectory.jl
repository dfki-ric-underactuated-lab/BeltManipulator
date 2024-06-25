"""
Simple trajectory type that contains states, controls and constraint forces, if index-3 DAE model is concerned
"""
mutable struct Trajectory
    
    const times::AbstractArray                  # time stamps
    const N::Int                                # number of knot points

    X::Vector{Vector{Float64}}                  # state variables
    U::Vector{Vector{Float64}}                  # controls


    function Trajectory(m::Int, n::Int, times::AbstractArray)
        
        N = length(times)

        X = [zeros(m) for _ = 1:N]
        # initialising controls with random values for iLQR
        U = [zeros(n) for _ = 1:N-1]
    
        return new(times, N, X, U)
    end
end

function Trajectory(mod::Model, times::AbstractArray)
    Trajectory(mod.m, mod.n, times)
end

function Trajectory(mod::Model, h::Float64, tf::Float64)
    Trajectory(mod, 0:h:tf)
end
