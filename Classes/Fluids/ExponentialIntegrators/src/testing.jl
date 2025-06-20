fixed_methods = 
    (
     ETDRK4(autodiff = false),
     HochOst4(autodiff = false),
     SBDF4(autodiff = false),
    )
fixed_names = 
    (
     "ETDRK4",
     "HochOst4",
     "SBDF4",
    )
adaptive_methods = 
    (
     RK4(),
     Tsit5(),
     KenCarp4(autodiff = false),
     KenCarp58(autodiff = false),
    )
adaptive_names = 
    (
     "RK4",
     "Tsit5",
     "KenCarp4",
     "KenCarp58",
    )

methods = (ETDRK4(autodiff = false),
           HochOst4(autodiff = false),
           RK4(),
           Tsit5(),
           SBDF4(autodiff = false),
           KenCarp3(autodiff = false),
           KenCarp4(autodiff = false),
           KenCarp58(autodiff = false),
          )

names = ("ETDRK4",
         "HochOst4",
         "RK4",
         "Tsit5",
         "SBDF4",
         "KenCarp3",
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

function errors_and_times_tolerance(problem, tru, tols; N = 5)
    errs = Vector{Float64}()
    times = Vector{Float64}()

    sizehint!(errs, length(tols))
    sizehint!(times, length(tols))

    for tol ∈ tols
        time = sum([(@elapsed problem(tol)) for i ∈ 1 : N]) / N
        push!(times, time)
        val = problem(tol)
        push!(errs, relative_error(val, tru))
    end

    errs[isnan.(errs) .|| isinf.(errs) .|| (errs .> 1)] .= 1
    errs[iszero.(errs)] .= 1e-16

    return (times, errs)
end

function errors_and_times_step_size(problem, tru, step_sizes; N = 5)
    errs = Vector{Float64}()
    times = Vector{Float64}()

    sizehint!(errs, length(step_sizes))
    sizehint!(times, length(step_sizes))

    for dt ∈ step_sizes
        time = sum([(@elapsed problem(dt)) for i ∈ 1 : N]) / N
        push!(times, time)
        val = problem(dt)
        push!(errs, relative_error(val, tru))
    end

    errs[isnan.(errs) .|| isinf.(errs) .|| (errs .> 1)] .= 1
    errs[iszero.(errs)] .= 1e-16

    return (times, errs)
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

function convergence_time_plot(times, errs, names)
    set_theme!(theme_latexfonts())
    fig = Figure()
    ax = Axis(fig[1,1], 
              xscale = log10, yscale = log10,
              xlabel = "Time (seconds)",
              ylabel = "Relative Error",
              xminorticksvisible = true,
              xminorgridvisible = true,
              xminorticks = IntervalsBetween(5)
             )
    for (time, err, name) ∈ zip(times, errs, names)
        lines!(ax, time, err, label = name)
    end

    Legend(fig[1,2], ax, "Method")

    colsize!(fig.layout, 1, Aspect(1,1))
    resize_to_layout!(fig)

    return fig
end

function error_plot(ivp, N_spatial, tspan, step_sizes, tols, tru_dt = 1e-5; N = 5, tru_solver = ETDRK4(autodiff = false))
    print("Getting true solution...")
    tru = ivp(N_spatial, tspan, tru_solver; save_everystep = false, dt = tru_dt, abstol = 1e-16, reltol = 1e-16)[3][end]
    println(" Done")
    
    errs  = Vector{Vector{Float64}}()
    times = Vector{Vector{Float64}}()
    println("\nRunning fixed time step tests:")
    for (method, name) ∈ zip(fixed_methods, fixed_names)
        println("\t" * name)
        problem(dt) = ivp(N_spatial, tspan, method; save_everystep = false, dt = dt)[3][end]
        time, err = errors_and_times_step_size(problem, tru, step_sizes; N = N)
        push!(times, time)
        push!(errs,  err)
    end
    println("Done")

    println("\nRunning adaptive time step tests:")
    for (method, name) ∈ zip(adaptive_methods, adaptive_names)
        println("\t" * name)
        problem(tol) = ivp(N_spatial, tspan, method; save_everystep = false, abstol = tol, reltol = tol)[3][end]
        time, err = errors_and_times_tolerance(problem, tru, tols, N = N)
        push!(times, time)
        push!(errs,  err)
    end
    println("Done")

    fig = convergence_time_plot(times, errs, (fixed_names..., adaptive_names...))

    return fig
end

function burgers_plot()
    tspan = (0, 1)
    N_tests    = 100 
    N_spatial  = 256
    step_sizes = 10 .^ range(-4, 0, N_tests)
    tols = 10 .^ range(-12, -1, N_tests)
    N_avg = 10

    fig = error_plot(burgers, N_spatial, tspan, step_sizes, tols; N = N_avg, tru_solver = Tsit5())
    fig.content[1].title = "Burgers Error"

    return fig
end

function schrodinger_plot()
    tspan = (0, 1)
    N_tests    = 50
    N_spatial  = 256
    step_sizes = 10 .^ range(-4, 0, N_tests)
    tols = 10 .^ range(-12, -1, N_tests)
    N_avg = 3

    fig = error_plot(schrodinger, N_spatial, tspan, step_sizes, tols; N = N_avg)
    fig.content[1].title = "Schrodinger's Equation Error"

    return fig
end

function kuramoto_plot()
    tspan = (0, 30)
    N_tests    = 50
    N_spatial  = 128
    step_sizes = 10 .^ range(-3, 0, N_tests)
    tols = 10 .^ range(-12, -1, N_tests)
    N_avg = 3

    fig = error_plot(kuramoto, N_spatial, tspan, step_sizes, tols; N = N_avg)
    fig.content[1].title = "Kuramoto Sivashinsky Equation Error"

    return fig
end

function KdV_plot()
    tspan = (0, 1e-3)
    N_tests    = 50
    N_spatial  = 512
    step_sizes = 10 .^ range(-5.5, -4, N_tests)
    tols = 10 .^ range(-8, -1, N_tests)
    N_avg = 3

    fig = error_plot(KdV, N_spatial, tspan, step_sizes, tols, 1e-6; N = N_avg)
    fig.content[1].title = "KdV Equation Error"

    return fig
end
