#=
iLQR trajectory for strawberry picking
=#
using BeltManipulator

begin
    urdf_path = abspath(joinpath(@__DIR__, "../", "urdf/shivaa_manipulator_local_copy.urdf"))

    # teeth numbers of gear stages
    n1 = 16; n2 = 16; n3 = 24; n4 = 24; n5i = 48; n5o = 24; n6 = 48; n7i = 48; n7o = 28; n8 = 42; n9i = 24; n9o = 20; n10i = 24; n10o = 20; n11 = 20; n12 = 20;

    # if no gears are considered
    # n1 = 1; n2 = 1; n3 = 1; n4 = 1; n5i = 1; n5o = 1; n6 = 1; n7i = 1; n7o = 1; n8 = 1; n9i = 1; n9o = 1; n10i = 1; n10o = 1; n11 = 1; n12 = 1;

    # building structure matrix
    g1 = (n1/n5i)*(n5o/n6)
    g2 = (n2/n7i)*(n7o/n8)
    g3 = n3/n9i
    g4 = n4/n10i
    g5 = (n3/n9i)*(n9o/n11)
    g6 = (n4/n10i)*(n10o/n12)

    S = [[1/g1, 1, 1, 1] [0, 1/g2, 1/g3, 1/g4] [0, 0, 1/g5, 1/g6] [0, 0, 1/g5, -1/g6]]

    # friction is not considered for the trajectory comparison
    # friction_coeffs = [(1.1, 0.1), (0.5, 0.05), (0.2, 0.02), (0.2, 0.02)]

    arm_js = RBDMechanism(urdf_path, S = diagm(ones(4)), visualizer = false);
    arm_as = RBDMechanism(urdf_path, friction = nothing, S = S, visualizer = false);

    # building the two four-DOF models
    model_js = ODEModel(4, arm_js);
    model_as = ODEModel(4, arm_as);
end 

# time-step and final time
h = 0.001
tf = 0.8

# start- and end-position
q0 = [0, 1.2π/2, 0, 0]
qf = [-π, -π/2, -π/2, 0]
 
# checking poses
# MeshCatMechanisms.set_configuration!(arm_js.mvis, qf)

### In joint-space
begin
    trj_joptimal_js = Trajectory(model_js, h, tf);

    op_js = iLQRParameters(
        x0 = [q0; [0, 0, 0, 0]], 
        xf = [qf; [0, 0, 0, 0]], 
        model = model_js, 
        times = trj_joptimal_js.times, 
        Q = diagm([1e-1*ones(4); 1e-1*ones(4)]), 
        R = diagm([1e-1, 1e-1, 1e-1, 5e2]), 
        Qf = diagm([[6e4, 3e4, 2e4, 1e3]; 1e4ones(4)]), 
        tol = 1e-3,
        reg_limit = 50,
        dβ = 1000.,
        a_red = 0.5,
        regularization = :state,
        second_order = false,
        i_max = 200
    );

    iLQR!(trj_joptimal_js, model_js, op_js)
end

### In actuation-space
begin
    trj_aoptimal_as = Trajectory(model_as, h, tf);

    op_as = iLQRParameters(
        x0 = [q0; [0, 0, 0, 0]], 
        xf = [qf; [0, 0, 0, 0]], 
        model = model_as, 
        times = trj_aoptimal_as.times,  
        ### weights for arm with gear ratios
        Q = diagm([1e-1*ones(4); 1e-1*ones(4)]), 
        R = diagm([1e1, 1e1, 1e1, 1e1]), 
        Qf = diagm([[3e5, 4e4, 3e4, 1e3]; 1e4ones(4)]),
        ### 
        ### weights for arm without gear ratios
        # Q = diagm([1e-1*ones(4); 1e-1*ones(4)]), 
        # R = diagm([1e-1, 1e-1, 1e1, 1e1]), 
        # Qf = diagm([[6e4, 1e4, 3e4, 1e3]; 1e4ones(4)]), 
        ###
        tol = 1e-3,
        reg_limit = 50,
        dβ = 1000.,
        a_red = 0.5,
        regularization = :state,
        second_order = false,
        i_max = 200
    );

    iLQR!(trj_aoptimal_as, model_as, op_as)
end


# mapping states of the actuation-space optimal trajectory in actuation space
joint2act_map!(trj_aoptimal_as, arm_as, states = true, torque = false);

# creating naive trajectories in both spaces
trj_naive_js = Trajectory(model_js, h, tf);
naive_trajectory!(trj_naive_js, arm_js, op_js);
trj_naive_as = deepcopy(trj_naive_js);
joint2act_map!(trj_naive_as, arm_as, states = true, torque = true);

# creating comparative trajectories and map them into the respective spaces
trj_joptimal_as = deepcopy(trj_joptimal_js);
joint2act_map!(trj_joptimal_as, arm_as, states = true, torque = true);
trj_aoptimal_js = deepcopy(trj_aoptimal_as);
act2joint_map!(trj_aoptimal_js, arm_as, states = true, torque = true);

# comparing everything
e_fig = BeltManipulator.power_comparison([
    (trj_naive_js, trj_naive_as),
    (trj_joptimal_js, trj_joptimal_as),
    (trj_aoptimal_js, trj_aoptimal_as)], interactive = false)

BeltManipulator.CairoMakie.save("power_comparison.png", e_fig)


# # Code for creating the overlayed graphics
# # before usage: SET DOWNLOAD DIRECTORY OF BROWSER TO 'meshcat_plots' 
# rate = 20
# n = length(1:rate:length(trj_naive_js.X)) 
# vis_frame_saver(arm_as, trj_naive_js.X, 1:rate:length(trj_naive_js.X))
# plt = frame_overlay(pwd()*"/images/meshcat_plots/",n,alpha_start=0.5) 
# savefig(pwd()*"/images/naive_traj.png")
# display(plt)  
