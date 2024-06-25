"""
    commutation_matrix(m::Int, n::Int)
    
Creates a commutation matrix as SparseMatrixCSC, like described in https://en.wikipedia.org/wiki/Commutation_matrix
"""
function commutation_matrix(m::Int, n::Int)
    sparse(1:m*n, [i + j*m for i = 1:m for j = 0:n-1], ones(m*n))
end

"""
Copies the state from Mechanism into the k-th entry of the Trajectory type
"""
function store_results!(trj::Trajectory, mec::Mechanism, k::Int)
    trj.X[k][1:mec.m] .= mec.q
    trj.X[k][mec.m+1:end] .= mec.q̇
end

"""
Checks, whether matrix can be inverted by simple determinant
"""
isinvertible(A::Matrix{Float64}) = !isapprox(det((A)), 0, atol = 1e-8)

"""
    function chol!(A::AbstractMatrix)

Splitting stuff up into lightweight posdef check
"""
function chol!(A::AbstractMatrix)
    C, info = LinearAlgebra._chol!(A, LowerTriangular)
    Cholesky(C.data, 'L', info)
end

"""
Returns moving average of array vs. Mean value is calucated using n points  
"""
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

"""
Symmetric logarithmic function - see J. B. W. Webber, “A bi-symmetric log 
transformation for wide-range data,” vol. 24, p. 027001, dec 2012.
"""
function symmetric_log(x::Real, c_factor::Real = 1/log(10))
    return sign(x)*log(1 + abs(x/c_factor))
end

"""
    copy_traj!(trjA::Trajectory, trjB::Trajectory)

Copies the states and controls of trajectory B to A and by that assumes that the have then same length and temporal distances
"""
function copy_traj!(trjA::Trajectory, trjB::Trajectory)
    copyto!(trjA.X, trjB.X)
    copyto!(trjA.U, trjB.U)
end

"""
Computes numerically the derivative of a state array by the underlying time-stamps
"""
function numerical_time_derivative(state::AbstractArray, times::AbstractArray)
    @assert length(state) == length(times) "Lengths of arrays does not match"

    state_diff = deepcopy(state)
    for i = 1:length(state)-1
        state_diff[i] = (state[i+1] - state[i])/(times[i+1] - times[i])
    end
    return state_diff
end

"""
    trapezoidal_rule(y::AbstractArray, x::AbstractArray)

Gives the approximated integral for the discrete arrays 'x' and 'y' by the trapezoidal rule.
"""
function trapezoidal_rule(y::AbstractArray, x::AbstractArray)
    @assert length(x) == length(y) "x and y must be of same length!"
    return sum([(y[k-1] + y[k])/2*(x[k] - x[k-1]) for k = 2:length(x)])
end

"""
    create_local_urdf_copy(urdf::String)

Creates a copy of a URDF with the correct relative paths to the mesh files in the package and assumes that 'urdf' is located in the path 'urdf/'
"""
function create_local_urdf_copy(urdf::String)
    # altering the relative paths in URDF-file, once this package is loaded
    filename = abspath(joinpath(@__DIR__, "..", "urdf/$(urdf).urdf"))
    xdoc = parse_file(filename)
    xroot = root(xdoc)
    # changing attributes only of the visual - there is also collision!
    for link in get_elements_by_tagname(xroot, "link")
        mesh = get_elements_by_tagname(get_elements_by_tagname(get_elements_by_tagname(link, "visual")[1], "geometry")[1], "mesh")[1]
        value = attribute(mesh, "filename")[2:end]
        ### NEVER USE JOINPATH HERE!!!!
        new_value = abspath(String(@__DIR__)*"/../"*value)
        set_attribute(mesh, "filename", new_value)
    end
    # save file as copy in package
    save_file(xdoc, joinpath(@__DIR__, "..", "urdf/$(urdf)_local_copy.urdf"))
end  


"""
    joint2act_map!(trj::Trajectory, mec::RBDMechanism; states::Bool = true, torque::Bool = true)

Translates a Trajectory in joint-space into actuation-space
"""
function joint2act_map!(trj::Trajectory, mec::RBDMechanism; states::Bool = true, torque::Bool = true)
    if states == torque == false
        @warn "Nothing is mapped"
    end

    if states
        for i = 1:trj.N
            trj.X[i] = [mec.S zeros(size(mec.S)); zeros(size(mec.S)) mec.S]*trj.X[i]
        end
    end
    if torque
        for i = 1:trj.N-1
            trj.U[i] = mec.S_'*trj.U[i]
        end
    end
end

"""
    act2joint_map!(trj::Trajectory, mec::RBDMechanism; states::Bool = true, torque::Bool = true)

Translates a Trajectory in actuation-space into joint-space
"""
function act2joint_map!(trj::Trajectory, mec::RBDMechanism; states::Bool = true, torque::Bool = true)
    if states == torque == false
        @warn "Nothing is mapped"
    end

    if states
        for i = 1:trj.N
            trj.X[i] = [mec.S_ zeros(size(mec.S)); zeros(size(mec.S)) mec.S_]*trj.X[i]
        end
    end
    if torque
        for i = 1:trj.N-1
            trj.U[i] = mec.S'*trj.U[i]
        end
    end
end

"""
    append(trj1::Trajectory, trj2::Trajectory)

Returns two appended trajectories by casting all their fields together
"""
function append(trj1::Trajectory, trj2::Trajectory)
    return Trajectory([Vector(trj1.times); Vector(trj2.times[2:end]) .+ trj1.times[end]], [trj1.X; trj2.X[2:end]], [trj1.U; trj2.U])
end