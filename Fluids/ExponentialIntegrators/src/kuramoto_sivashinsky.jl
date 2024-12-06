function kuramoto(N = 129, tspan = (0, 150), method = ETDRK4(autodiff = false); kwargs... )
    order = N ÷ 2                                   # Max complex exponential order
    k = [0:order; -(order - 1 + N % 2):-1] / 16     # Wave numbers

    xlims = (0, 32π * (1 - 1/N))                    # Ignore endpoint
    x = range(xlims..., N)                          # Grid points

    u0 = cos.(x / 16) .* (1 .+ sin.(x / 16))        # Initial condition
    Fu0 = fft(u0)                                   # Frequency domain initial condition
     
    L = DiagonalOperator(k.^2 - k.^4)               # Linear part
    Non(u, p, t) = -0.5im * k .* fft(real(ifft(u)).^2) # Nonlinear part
    
    problem = SplitODEProblem(L, Non, Fu0, tspan)   # Define split u' = Lu + Nu problem
    sol = solve(problem, method; kwargs...)         # Solve the system

    u = real.(ifft.(sol.u))                         # Convert solution back to real space

    return (sol.t, x, u)
end

function plot_kuramoto(t, x, u)
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

function animate_kuramoto(t, x, u; filename = "animation")
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

# function kuramoto(;N = 129, tspan = (0, 65), method = ETDRK4, Δt = .1)
#     order = N ÷ 2                                   # Max complex exponential order
#     k = (-(order):(order - (1 - (N % 2)))) / 16     # Complex exponential orders
#     # k = [(-(order):(order - (1 - (N % 2)))) / 16...]     # Complex exponential orders
#     # k[1] = 0
     
#     L = DiagonalOperator(k.^2 - k.^4)               # Linear part

#     # Non(u, p, t) = -convolve(u, im * k .* u)[1 : N] # Nonlinear part
#     Non(u, p, t) = -im * k .* fftshift(fft(real(ifft(ifftshift(u))).^2))
    
#     # Initial conditions
#     xlims = (0, 32π * (1 - 1/N))                    # Ignore endpoint for interpolation
#     # xlims = (32π/N, 32π) 
#     x = range(xlims..., N)                          # Grid points
#     u0 = cos.(x / 16) .* (1 .+ sin.(x / 16))        # Initial condition
#     # Fu0 = fftshift(fft(u0)) / N                     # Frequency domain initial condition
#     Fu0 = fftshift(fft(u0))                         # Frequency domain initial condition

#     problem = SplitODEProblem(L, Non, Fu0, tspan)   # Define split u' = Lu + Nu problem
#     sol = solve(problem,                            # Solve system
#                 method(autodiff = false), 
#                 dt = Δt)

#     # U = hcat(real.(N * ifft.(ifftshift.(sol.u)))...)
#     U = hcat(real.(ifft.(ifftshift.(sol.u)))...)
#     return (sol.t, x, U)
# end
