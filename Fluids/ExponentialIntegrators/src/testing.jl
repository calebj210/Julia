methods = (ETDRK4(autodiff = false),
           HochOst4(autodiff = false),
           RK4(),
           Tsit5(),
           SBDF4(autodiff = false),
           KenCarp4(autodiff = false),
           KenCarp58(autodiff = false),
          )

names = ("ETDRK4",
         "HochOst4",
         "RK4",
         "Tsit5",
         "SBDF4",
         "KenCarp4",
         "KenCarp58",
        )

function relative_error(val, tru; N = Inf)
    return norm(val - tru, N) / norm(tru, N)
end

function errors_and_times(problem, tru, step_sizes)
    errs = Vector{Float64}()
    times = Vector{Float64}()
    sizehint!(errs, length(step_sizes))
    sizehint!(times, length(step_sizes))

    for dt ∈ step_sizes
        push!(times, @elapsed val = problem(dt))
        push!(errs, relative_error(val, tru))
    end

    errs[isnan.(errs) .|| isinf.(errs) .|| (errs .> 1)] .= 1
    errs[iszero.(errs)] .= 1e-16

    return (errs, times)
end

function convergence_plot(step_sizes, errs, times, names)
    set_theme!(theme_latexfonts())
    fig = Figure()
    ax1 = Axis(fig[1,1], 
              xreversed = true, 
              xscale = log10, yscale = log10,
              xlabel = L"Step Size $\Delta t$",
              ylabel = "Relative Error",
              xminorticksvisible = true,
              xminorgridvisible = true,
              xminorticks = IntervalsBetween(5)
             )
    ax2 = Axis(fig[1,2], 
              xreversed = true, 
              xscale = log10, yscale = log10,
              xlabel = L"Step Size $\Delta t$",
              ylabel = "Solve Time (seconds)",
              title  = "Time of Computation",
              xminorticksvisible = true,
              xminorgridvisible = true,
              xminorticks = IntervalsBetween(5)
             )
    for (err, time, name) ∈ zip(errs, times, names)
        lines!(ax1, step_sizes, err, label = name)
        lines!(ax2, step_sizes, time, label = name)
    end

    Legend(fig[1,3], ax1, "Method")

    colsize!(fig.layout, 1, Aspect(1,1))
    colsize!(fig.layout, 2, Aspect(1,1))
    resize_to_layout!(fig)

    return fig
end

function burgers_plot()
    tspan = (0, 1)
    N_time     = 50
    N_spatial  = 256
    step_sizes = 10 .^ range(-4, 0, N_time)

    print("Getting true solution...")
    tru = burgers(N_spatial, tspan, ETDRK4(autodiff = false); save_everystep = false, dt = 1e-5)[3][end]
    println(" Done")
    
    errs  = Vector{Vector{Float64}}()
    times = Vector{Vector{Float64}}()
    println("\nRunning tests")
    for (method, name) ∈ zip(methods, names)
        println("\t" * name)
        problem(dt) = burgers(N_spatial, tspan, method; save_everystep = false, dt = dt)[3][end]
        err, time = errors_and_times(problem, tru, step_sizes)
        push!(errs,  err)
        push!(times, time)
    end
    println("Done")

    fig = convergence_plot(step_sizes, errs, times, names)
    fig.content[1].title = "Burgers Error"

    return fig
end

function schrodinger_plot()
    tspan = (0, 1)
    N_time     = 50
    N_spatial  = 256
    step_sizes = 10 .^ range(-4, 0, N_time)

    println("Getting true solution")
    tru = schrodinger(N_spatial, tspan, ETDRK4(autodiff = false); save_everystep = false, dt = 1e-5)[3][end]
    
    errs  = Vector{Vector{Float64}}()
    times = Vector{Vector{Float64}}()
    println("\nRunning tests")
    for (method, name) ∈ zip(methods, names)
        println("\t" * name)
        problem(dt) = schrodinger(N_spatial, tspan, method; save_everystep = false, dt = dt)[3][end]
        err, time = errors_and_times(problem, tru, step_sizes)
        push!(errs,  err)
        push!(times, time)
    end
    println("Done")

    fig = convergence_plot(step_sizes, errs, times, names)
    fig.content[1].title = "Schrodinger's Equation Error"

    return fig
end

function kuramoto_plot()
    tspan = (0, 30)
    N_time     = 50
    N_spatial  = 128
    step_sizes = 10 .^ range(-3, 0, N_time)

    println("Getting true solution")
    tru = kuramoto(N_spatial, tspan, ETDRK4(autodiff = false); save_everystep = false, dt = 1e-4)[3][end]
    
    errs  = Vector{Vector{Float64}}()
    times = Vector{Vector{Float64}}()
    println("\nRunning tests")
    for (method, name) ∈ zip(methods, names)
        println("\t" * name)
        problem(dt) = kuramoto(N_spatial, tspan, method; save_everystep = false, dt = dt)[3][end]
        err, time = errors_and_times(problem, tru, step_sizes)
        push!(errs,  err)
        push!(times, time)
    end
    println("Done")

    fig = convergence_plot(step_sizes, errs, times, names)
    fig.content[1].title = "Kuramoto Sivashinsky Equation Error"

    return fig
end
