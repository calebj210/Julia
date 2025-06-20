function burgers(N = 128, tspan = (0, 1), method = ETDRK4(autodiff = false); kwargs...)
    λ = 0.03

    order = N ÷ 2                                   # Max complex exponential order
    k = [0:order; -(order - 1 + N % 2):-1]          # Wave numbers

    xlims = (-π, π * (1 - 2/N))                     # Ignore endpoint
    x = range(xlims..., N)                          # Grid points

    u0 = exp.(-10sin.(x / 2).^2)                    # Initial condition
    Fu0 = fft(u0)                                   # Frequency domain initial condition
     
    L = DiagonalOperator(-λ * k.^2)                  # Linear part
    NL(u, p, t) = -.5im * k .* fft(real(ifft(u)).^2) # Nonlinear part
    
    problem = SplitODEProblem(L, NL, Fu0, tspan)    # Define split u' = Lu + Nu problem
    sol = solve(problem, method; kwargs...)         # Solve the system

    u = real.(ifft.(sol.u))                         # Convert solution back to real space

    return (sol.t, x, u)
end

function burgers_heatmap(N = 128, tspan = (0, 1), method = ETDRK4(autodiff = false); kwargs...)
    t, x, u = burgers(N, tspan, method; kwargs...)

    U = hcat(u...)

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax = Axis(fig[1,1],
              xlabel = "x",
              ylabel = "t"
             )

    plt = heatmap!(ax, x, t, U)

    Colorbar(fig[1,2], plt)

    colsize!(fig.layout, 1, Aspect(1,1))
    resize_to_layout!(fig)

    return fig
end

function burgers_animation(filename = "burgers", N = 128, tspan = (0,1), method = ETDRK4(autodiff = false); kwargs...)
    t, x, u = burgers(N, tspan, method; kwargs...)

    frame = Observable(1)
    U = @lift(u[$frame])
    u_matrix = hcat(u...)

    xlims = (minimum(minimum.(u)), maximum(maximum.(u)))

    fig = Figure()
    ax1 = Axis(fig[1,1], 
              title = @lift("t = $(round(t[$frame], digits = 2))"), 
              limits = (nothing, xlims),
              xlabel = L"x",
              ylabel = L"u")
    lines!(ax1, x, U, color = :orange)

    ax2 = Axis(fig[1,2], xlabel = L"x", ylabel = "Time")
    plt = heatmap!(ax2, x, t, u_matrix, interpolate = true)
    Colorbar(fig[1,3], plt)
    lines!(ax2, x, @lift(repeat([t[$frame]], length(x))), color = :orange, linewidth = 2)

    colsize!(fig.layout, 1, Aspect(1,1))
    colsize!(fig.layout, 2, Aspect(1,1))
    resize_to_layout!(fig)

    framerate = 30
    timestamps = 1:length(t)

    record(fig, filename*".mp4", timestamps; framerate = framerate) do i
        frame[] = i
    end
end
