mutable struct iLQRParameters
    
    # general parameters
    const x0::Vector{Float64}                       # initial configuration
    const xf::Vector{Float64}                       # final configuration

    # optimal control parameters
    const Q::Matrix{Float64}
    const R::Matrix{Float64}
    const Qf::Matrix{Float64}

    const tol::Float64                              # stopping criteria for biggest term in feed-forward vector
    const i_max::Int                                # maximum overall iterations 
    const ls_max::Int                               # maximum line search trials
    const fail_max::Int                             # maximum number of subsequent fails allowed
    const second_order::Bool                        # controls the evaluation of second order derivatives
    regularization::Symbol                          # defines if and in which space regularization takes place
    z_low::Float64                                  # lower bound of expected cost reduction
    z_upp::Float64                                  # upper bound of expected cost reduction
    a_red::Float64                                  # reduction coefficient in line search

    β::Float64                                      # regularization
    const dβ::Float64                               # regularization increment
    const reg_max::Int                              # maximum number of repeated regularizations in one pass

    # new state and control values upon iteration
    trj_new::Trajectory

    # gains and feed-forward values for Trajectory
    K::Vector{Matrix{Float64}}
    d::Vector{Vector{Float64}}

    count::Int                              # counting variable for iterations
    fails::Int                              # failing variable

    """
        iLQRParameters(x0::Vector, xf::Vector, times::AbstractArray; Q::Matrix, R::Matrix, Qf::Matrix, tol::Float64, i_max::Int = 50)
    
    Standard constructor with predefined final time and given time array
    """
    function iLQRParameters(;x0::Vector, xf::Vector, model::Model, times::AbstractArray, Q::Matrix, R::Matrix, Qf::Matrix, tol::Float64, i_max::Int = 100, fail_max::Int = 5, dβ::Float64 = 10., reg_limit::Int = 10, ls_max::Int = 20, second_order::Bool = false, regularization::Symbol = :state, z_low::Float64 = 1e-4, z_upp::Float64 = 10., a_red::Float64 = 0.5)

        if size(Qf) != size(Q)
            throw(error("Final and intermediate cost matrix dimensions do not match. Qf ∈ ℝ^{$(size(Qf)[1]) × $(size(Qf)[2])} but Q ∈ ℝ^{$(size(Q)[1]) × $(size(Q)[2])}"))
        elseif size(Q)[1] != model.m && size(R)[1] != model.n
            throw(error("Mismatch in cost-matrix dimension and model dimension"))
        end 

        @assert regularization in([:state, :control]) "Keyword for regularization can only be ':state', or ':control'"

        @assert a_red < 1. "Reduction value needs to be < 1"

        # trajectory used for line search
        trj_new = Trajectory(model, times)    
        
        K = Vector{Matrix{Float64}}(fill(zeros(size(R)[1], size(Q)[1]), trj_new.N-1))
        d = Vector{Vector{Float64}}(fill(zeros(size(R)[1]), trj_new.N-1))

        new(x0, xf, Q, R, Qf, tol, i_max, ls_max, fail_max, second_order, regularization, z_low, z_upp, a_red, 0., dβ, reg_limit, trj_new, K, d, 0, 0)
    end

end