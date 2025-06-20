function KdV(N = 512, tspan = (0, .1), method = ETDRK4(autodiff = false); kwargs... )
    order = N ÷ 2                                   # Max complex exponential order
    k = [0:order; -(order - 1 + N % 2):-1]          # Wave numbers
    # k = [0:order; -(order - 1 + N % 2):-1] * π      # Wave numbers

    xlims = (-π, π * (1 - 2/N))                     # Ignore endpoint
    # xlims = (0, 2 * (1 - 1/N))                     # Ignore endpoint
    x = range(xlims..., N)                          # Grid points

    A = 25; B = 16
    u0 = 3 * A^2 * sech.(A * (x .+ 2) / 2).^2 +
         3 * B^2 * sech.(B * (x .+ 1) / 2).^2       # Initial condition
    # u0 = cos.(π * x)
    Fu0 = fft(u0)                                   # Frequency domain initial condition
     
    L = DiagonalOperator(im * k.^3)                 # Linear part
    # L = DiagonalOperator(0.022^2 * im * k.^3)       # Linear part
    Non(u, p, t) = -0.5im * k .* fft(real(ifft(u)).^2) # Nonlinear part
    
    problem = SplitODEProblem(L, Non, Fu0, tspan)   # Define split u' = Lu + Nu problem
    sol = solve(problem, method; kwargs...)         # Solve the system

    u = real.(ifft.(sol.u))                         # Convert solution back to real space

    return (sol.t, x, u)
end

function plot_KdV(t, x, u)
    T = range(t[1], t[end], length(t))
    X = x
    U = hcat(u...)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel = L"x", ylabel = "Time")
    plt = heatmap!(ax, X, T, U, interpolate = true)
    Colorbar(fig[1,2], plt)
    colsize!(fig.layout, 1, Aspect(1,1))
    resize_to_layout!(fig)

    return fig
end

function animate_KdV(t, x, u; filename = "animation")
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

    framerate = 60
    timestamps = 1:length(t)

    record(fig, filename*".mp4", timestamps; framerate = framerate) do i
        frame[] = i
    end
end
