module BeltManipulator


using LinearAlgebra
using SparseArrays
export diagm

using ForwardDiff

using Dates
using FileIO
import LightXML: parse_file, root, get_elements_by_tagname, attribute, set_attribute, save_file

using CairoMakie
# this is used to selectively display interactive plots with OpenGL
import GLMakie as GL
# Ensure standard backend Cairo is activated  
CairoMakie.activate!()  

using MeshCat
using MeshCatMechanisms
using RigidBodyDynamics

const FD = ForwardDiff

# basic types for the creation of mechanisms and specifically, parllel mechanisms
include("./Types/RBDMechanism.jl")
export RBDMechanism
include("./Types/Model.jl")
export ODEModel
include("./Types/Trajectory.jl")
export Trajectory

# hyperparameters for trajectory optimization and simulation
include("./Types/iLQRParameters.jl")
export iLQRParameters
include("./Types/ControlParameters.jl")
export ControlParameters

include("./discrete_dynamics.jl")
export runge_kutta4, naive_trajectory!, torques_from_trajectory!

include("./mechanism_dynamics.jl")

include("./controls.jl")

include("./visualization.jl")
export show_config, animate_mechanism!, plot_trajectory, power_comparison, vis_frame_saver, frame_overlay

include("./ilqr.jl")
export iLQR!

# package internal helper functions
include("./utils.jl")
export symmetric_log, joint2act_map!, act2joint_map!, append, create_local_urdf_copy

# this creates a copy with the correct paths to the mesh files
create_local_urdf_copy("shivaa_manipulator")

end # module BeltManipulator
