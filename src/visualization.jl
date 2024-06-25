"""
    visual_backend_selector(bend::String)

Allows for selecting the backend 'bend' for the graphics package makie.
This is useful for temporary allowing interactive plotting.       
"""
function visual_backend_selector(bend::String)
    if bend=="Cairo"
        CairoMakie.activate!()
    elseif bend=="OpenGL"
        GL.activate!()
    else 
        throw(error("Valid backends are 'Cairo' or 'OpenGL'"))
    end 
end  

"""
MeshCat animation of the manipulator
"""
function animate_mechanism!(mec::RBDMechanism, trj::Trajectory; step::Int = 1)
    # only taking positions from state array
    xs = [trj.X[i][1:num_positions(mec.rbd_model)] for i = 1:step:length(trj.times)]
    animation = MeshCat.Animation(mec.mvis, trj.times[begin:step:end], xs);
    setanimation!(mec.mvis, animation);
end

"""
    plot_trajectory(trj::Trajectory; cscheme=:Dark2_8, size = (1200, 400), acceleration::Bool = false, title = "", interactive = true)
    
Plots states and controls of a given trajectory
"""
function plot_trajectory(trj::Trajectory; cscheme=:Dark2_8, size = (1200, 400), acceleration::Bool = false, title = "", interactive = true)
    if interactive 
        visual_backend_selector("OpenGL") 
    else 
        visual_backend_selector("Cairo") 
    end 
    fig = Figure(fontsize = 24)
    ax1 = Axis(fig[1,1],ylabel = "x, v [rad, rad/sec]",width=size[1],height=size[2], title = title,xlabel = "time [s]")

    n = Int(length(trj.X[1])/2) # half number of states
    
    plots_pos = []
    plots_vel = []
    plots_acc = []

    for i = 1:n
        ppos = lines!(ax1, trj.times, [trj.X[j][i] for j = 1:trj.N], linewidth = 2, color = i, colormap = cscheme, colorrange = (1, 8), label = "pos $i")
        pvel = lines!(ax1, trj.times, [trj.X[j][n+i] for j = 1:trj.N], linewidth = 2, linestyle = :dash, color = i, colormap = cscheme, colorrange = (1, 8), label = "vel $i")
        push!(plots_pos,ppos)
        push!(plots_vel,pvel)        
        if acceleration
            pacc = lines!(trj.times, numerical_time_derivative([trj.X[j][n+i] for j = 1:trj.N], trj.times), linewidth = 2, linestyle = :dashdot, color = i, colormap = cscheme, colorrange = (1, 8), label = "acc $i")
            push!(plots_acc, pacc)
        end
    end

    fig[1,2] = Legend(fig,[plots_pos;plots_vel;plots_acc], [["pos $i" for i=1:n];["vel $i" for i=1:n];["acc $i" for i=1:length(plots_acc)]], "state trajectory")

    ax2 = Axis(fig[2,1],xlabel = "time [s]", ylabel = "control [Nm]", width=size[1],height=size[2])
    m = length(trj.U[1]) # number of controls
    plots_ctrl = []
    for i = 1:m
        pctrl = lines!(ax2, trj.times[1:end-1], [trj.U[j][i] for j = 1:trj.N-1], linewidth = 3, alpha = 0.6, label = "mot $i", color = i, colormap = cscheme, colorrange = (1, 8),)
        push!(plots_ctrl,pctrl)
    end
    fig[2,2] = Legend(fig,plots_ctrl, ["mot $i" for i=1:m], "input trajectory")
    if interactive DataInspector(fig) end 
    resize_to_layout!(fig)
    fig         
    return fig
end 

"""
    plot_policy(op::ILQRParameters, cscheme=:Dark2_8, size = (1600, 400), title = ""; interactive = true)

Shows the feed-forward and feed-back values of a current optimization stage
"""
function plot_policy(op::iLQRParameters, cscheme=:Dark2_8, size = (1600, 400), title = ""; interactive = true)
    if interactive 
        visual_backend_selector("OpenGL") 
    else 
        visual_backend_selector("Cairo") 
    end 

    fig = Figure(fontsize = 24)
    ax1 = Axis(fig[1,1],ylabel = "p",width=size[1],height=size[2], title = title,xlabel = "time [s]")
    plots_ff = []    
    for i = 1:length(op.d[1])
        ff = lines!(ax1, op.trj_new.times[1:end-1], [op.d[j][i] for j = 1:op.trj_new.N-1], linewidth = 2, color = i, colormap = cscheme, colorrange = (1, 8), label = "ff $i")
        push!(plots_ff,ff)
    end
    fig[1,2] = Legend(fig,plots_ff, ["$i" for i=1:length(op.d[1])], "feed forward")

    ax2 = Axis(fig[2,1],ylabel = "σ(K)",width=size[1],height=size[2], title = title,xlabel = "time [s]")
    plots_fb=[]
    for i = 1:length(vec(op.K[1]))
        fb = lines!(ax2, op.trj_new.times[1:end-1], [vec(op.K[j])[i] for j = 1:op.trj_new.N-1], lw = 1, label = "vec(K) $i") # [svd(op.K[j]).S[i] for j = 1:op.N-1] when using svd
        push!(plots_fb,fb)
    end
    fig[2,2] = Legend(fig,plots_fb, ["$i" for i=1:length(vec(op.K[1]))], "feed backward")
    if interactive DataInspector(fig) end 
    resize_to_layout!(fig)
    fig         
    return fig
end

"""
    compare_trajectories(trj::Vector{Trajectory}, titles::Union{Nothing, Vector{String}} = nothing; figsize = (600, 400), cscheme=:Dark2_8, interactive::Bool = true)
    
FOR NOW ONLY 2 TRAJECTORIES POSSIBLE. Comparison by state, controls and position error.
"""
function compare_trajectories(trj::Vector{Trajectory}, titles::Union{Nothing, Vector{String}} = nothing; figsize = (600, 400), cscheme=:Dark2_8, interactive::Bool = true)
    if interactive 
        visual_backend_selector("OpenGL") 
    else 
        visual_backend_selector("Cairo") 
    end 


    fig = Figure(fontsize = 24)
    ax1 = Axis(fig[1,1], xlabel = "time [s]", ylabel = "x rad",width=figsize[1], height=figsize[2])
    ax2 = Axis(fig[2,1], xlabel = "time [s]", ylabel = "v rad/sec",width=figsize[1], height=figsize[2])
    ax3 = Axis(fig[1,3], xlabel = "time [s]", ylabel = "control [Nm]", width=figsize[1], height=figsize[2])
    ax4 = Axis(fig[2,3], xlabel = "time [s]", ylabel = "sym_log error", width=figsize[1], height=figsize[2])


    n = Int(length(trj[1].X[1])/2) # half number of states
    m = length(trj[1].U[1]) # number of controls
    
    plots_pos1 = []
    plots_pos2 = []

    plots_vel1 = []
    plots_vel2 = []

    plots_ctrl1 = []
    plots_ctrl2 = []

    plots_err = []                     

    for i = 1:n
        ppos1 = lines!(ax1, trj[1].times, [trj[1].X[k][i] for k = 1:trj[1].N], color = i, colormap = cscheme, colorrange = (1, 8), linewidth = 2,label = "pos $i")  #traj 1 
        ppos2 = lines!(ax1, trj[2].times, [trj[2].X[k][i] for k = 1:trj[2].N], color = i, colormap = cscheme, colorrange = (1, 8), linewidth = 2,label = "pos $i",linestyle = :dash)  #traj 2
        push!(plots_pos1, ppos1)
        push!(plots_pos2, ppos2)

        pvel1 = lines!(ax2, trj[1].times, [trj[1].X[k][n+i] for k = 1:trj[1].N], color = i, colormap = cscheme, colorrange = (1, 8), linewidth = 2, label = "vel $i")
        pvel2 = lines!(ax2, trj[2].times, [trj[2].X[k][n+i] for k = 1:trj[2].N], color = i, colormap = cscheme, colorrange = (1, 8), linewidth = 2, linestyle = :dash, label = "vel $i")
        push!(plots_vel1, pvel1)        
        push!(plots_vel2, pvel2)        

        perr = lines!(ax4, trj[2].times, [symmetric_log(trj[1].X[k][i] - trj[2].X[k][i]) for k = 1:trj[2].N], color = i, colormap = cscheme, colorrange = (1, 8), linewidth = 2, label = "err $i")   
        push!(plots_err, perr)
    end
    
    for i = 1:m
        pctrl1 = lines!(ax3, trj[1].times[1:end-1], [trj[1].U[k][i] for k = 1:trj[1].N-1], color = i, colormap = cscheme, colorrange = (1, 8), label = "mot $i")
        pctrl2 = lines!(ax3, trj[2].times[1:end-1], [trj[2].U[k][i] for k = 1:trj[2].N-1], color = i, colormap = cscheme, colorrange = (1, 8), label = "mot $i", linestyle = :dash)    
        push!(plots_ctrl1,pctrl1)
        push!(plots_ctrl2,pctrl2)
    end
    
    fig[1,2] = Legend(fig, [plots_pos1; plots_pos2], [["$(titles[1]) $i" for i=1:n];["$(titles[2]) $i" for i=1:n]], "position")
    fig[2,2] = Legend(fig, [plots_vel1; plots_vel2], [["$(titles[1]) $i" for i=1:n];["$(titles[2]) $i" for i=1:n]], "velocity")
    fig[1,4] = Legend(fig, [plots_ctrl1; plots_ctrl2], [["$(titles[1]) $i" for i=1:m];["$(titles[2]) $i" for i=1:m]], "input")
    fig[2,4] = Legend(fig, plots_err, ["$(titles[1]) $i" for i=1:n], "pos. error")

    if interactive DataInspector(fig) end 
    resize_to_layout!(fig)
    fig         
    return fig
end

"""
    plot_trajectory(trj::Trajectory; cscheme=:Dark2_8, size = (1200, 400), acceleration::Bool = false, title = "", interactive = true)
    
Power and energy losses of a trajectory
"""
function plot_power_use(trj::Trajectory; cscheme=:Dark2_8, size = (1200, 400), title = "", interactive = true)
    if interactive 
        visual_backend_selector("OpenGL") 
    else 
        visual_backend_selector("Cairo") 
    end 
    fig = Figure(fontsize = 24)
    ax1 = Axis(fig[1,1],ylabel = "power [W]",width=size[1],height=size[2], title = title,xlabel = "time [s]")
    ax2 = Axis(fig[2,1],xlabel = "time [s]", ylabel = "energy [J]", width=size[1],height=size[2])

    n = Int(length(trj.X[1])/2) # half number of states
    
    plots_pow = []
    plots_ene = []

    for i = 1:n
        # computing power in motors
        pow = [trj.X[k][n+i]*trj.U[k][i] for k = 1:trj.N-1]
        pow_l = [pow[k] < 0 ? abs(pow[k]) : 0. for k = 1:trj.N-1]
        ene = [trapezoidal_rule(pow[1:k], trj.times[1:k]) for k = 1:trj.N-1]
        ene_l = [trapezoidal_rule(pow_l[1:k], trj.times[1:k]) for k = 1:trj.N-1]
        powp = lines!(ax1, trj.times[1:end-1], pow, linewidth = 2, color = i, colormap = cscheme, colorrange = (1, 8), label = "mot $i")
        powlp = lines!(ax1, trj.times[1:end-1], pow_l, linewidth = 2, linestyle = :dash, color = i, colormap = cscheme, colorrange = (1, 8), label = "mot $i")
        enep = lines!(ax2, trj.times[1:end-1], ene, linewidth = 2, color = i, colormap = cscheme, colorrange = (1, 8), label = "mot $i")
        enelp = lines!(ax2, trj.times[1:end-1], ene_l, linewidth = 2, linestyle = :dash, color = i, colormap = cscheme, colorrange = (1, 8), label = "mot $i")
        push!(plots_pow, powp)
        push!(plots_pow, powlp)
        push!(plots_ene, enep)    
        push!(plots_ene, enelp)        
    end

    fig[1,2] = Legend(fig, plots_pow, [["pow $i" for i=1:n]; ["loss $i" for i=1:n]], "power")
    fig[2,2] = Legend(fig, plots_ene, [["ene $i" for i=1:n]; ["loss $i" for i=1:n]], "energy")
    
    if interactive DataInspector(fig) end 
    resize_to_layout!(fig)
    fig         
    return fig
end


"""
    power_comparison(trjs::Vector{Tuple{Trajectory, Trajectory}}; cscheme=:Dark2_8, size = (1200, 400), title = "", interactive = true) 
    
For the vector of Trajectory Tuples the power is compared. It is assumed that the first entry of the Tuple is a Trajectory in joint-space and the second in actuation-space.
"""
function power_comparison(trjs::Vector{Tuple{Trajectory, Trajectory}}; cscheme=:Dark2_8, size = (1200, 300), title = "", interactive = true)
    
    if interactive 
        visual_backend_selector("OpenGL") 
    else 
        visual_backend_selector("Cairo") 
    end 

    fig = Figure(fontsize = 18)
    line_axes = []
    power_axes = []
    loss_axes = []
    for i = 1:length(trjs)
        if i < length(trjs)
            push!(line_axes, Axis(fig[i,1], ylabel = "power [W]",width=size[1],height=size[2]))
        else
            push!(line_axes, Axis(fig[i,1], ylabel = "power [W]",width=size[1],height=size[2], xlabel = "time [s]"))
        end
        push!(power_axes, Axis(fig[i,2], xticks = (1:2, ["joint", "act."]), ylabel = "peak power [W]", width=size[1]/16, height=size[2]))
        push!(loss_axes, Axis(fig[i,3], xticks = (1:2, ["joint", "act."]), ylabel = "loss [J]", width=size[1]/16, height=size[2], yaxisposition = :right))
    end

    n = Int(length(trjs[1][1].X[1])/2) # half number of states
    space_colors = [RGBAf(255, 0, 0, 0.5), RGBAf(0.007, 1, 0.007, 0.5)]

    peak_power_scaling = []
    energy_loss_scaling = []

    # considering the different plots
    for (i, trj_) in enumerate(trjs)
        # taking joint- and actuation-space trajectories
        for (s, trj) in enumerate(trj_)
            max_powers = []
            energy_losses = []
            # looking at individual joints/motors
            for j = 1:n
                # computing power
                pow = [trj.X[k][n+j]*trj.U[k][j] for k = 1:trj.N-1]
                # computing power losses
                pow_l = [pow[k] < 0 ? -abs(pow[k]) : 0. for k = 1:trj.N-1]
                # computing energy loss at end of trajectory
                ene_l = trapezoidal_rule(pow_l, trj.times[1:end-1])
                # getting the maximal power of that joint/motor
                push!(max_powers, maximum(abs.(pow)))
                # add energy loss of that joint/motor
                push!(energy_losses, ene_l)

                # deciding label and appearance
                if i == 1 && s == 1
                    label_ = "joint $j"
                    linestyle_ = :dash
                elseif i == 1 && s == 2
                    label_ = "motor $j"
                    linestyle_ = :solid
                elseif i != 1 && s == 1
                    label_ = false
                    linestyle_ = :dash
                elseif i != 1 && s == 2
                    label_ = false
                    linestyle_ = :solid
                end
                
                # plotting things
                lines!(line_axes[i], trj.times[1:end-1], pow, linewidth = 2, color = j, colormap = cscheme, colorrange = (1, 8), linestyle = linestyle_, label = label_)
                band!(line_axes[i], trj.times[1:end-1], pow_l, zeros(length(ene_l)), color = (:red, 0.2))
                
            end
            # getting the maximal values out of the loop
            push!(peak_power_scaling, maximum(max_powers))
            push!(energy_loss_scaling, abs(sum(energy_losses)))

            # plotting peak power and losses
            barplot!(power_axes[i], [s], [maximum(max_powers)], color = space_colors[s], strokecolor = :black, strokewidth = 0.5)
            barplot!(loss_axes[i], [s], [abs(sum(energy_losses))], color = space_colors[s], strokecolor = :black, strokewidth = 0.5)
        end
    end
    # only putting legend on first plot
    axislegend(line_axes[1])
    
    # scaling the bar plots uniformly
    for i = 1:length(trjs)
        limits!(power_axes[i], nothing, nothing, 0, 1.05*maximum(peak_power_scaling))
        limits!(loss_axes[i], nothing, nothing, 0, 1.05*maximum(energy_loss_scaling))
    end

    if interactive DataInspector(fig) end 
    resize_to_layout!(fig)
    fig         
    return fig
end

"""
Saves screenshots of specified trajectory, adjust properties beforehand in the browser
"""
function vis_frame_saver(mech::RBDMechanism, X::AbstractArray, indices)
    if maximum(indices) > length(X)
        throw(error("indices are out of trajectory bounds"))
    end 
    for i in indices
        set_configuration!(mech.mvis, X[i][1:4])
        MeshCat.save_image(mech.vis)
    end 
end 

"""
Scales transparency by chaning alpha channel 
"""
scaletransparency(image, scale) = RGBA.(RGB.(image), scale.*alpha.(image))


"""
Create overlay plot with alpha gradient 
"""
function frame_overlay(path, n_files; flipdir=false, alpha_start=0, alpha_end=1)
    alpha = range(alpha_start,alpha_end,length=n_files)
    if flipdir
        alpha = reverse(alpha)
    end     
    img = FileIO.load(path*"meshcat.png")    
    img=scaletransparency(img, alpha[1])
    plt =  Plots.plot(img, size=(1200,800),xticks=false,yticks=false,axis=([], false))   
    for i in eachindex(alpha[2:end])
        img=FileIO.load(path*"meshcat($i).png")
        img=scaletransparency(img, alpha[i])
        Plots.plot!(img)        
    end 
    return plt
end  
