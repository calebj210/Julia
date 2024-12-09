function schrodinger(N = 128, tspan = (0, 1), method = ETDRK4(autodiff = false); kwargs...)
    order = N ÷ 2                                   # Max complex exponential order
    k = [0:order; -(order - 1 + N % 2):-1]          # Wave numbers

    xlims = (-π, π * (1 - 2/N))                     # Ignore endpoint
    x = range(xlims..., N)                          # Grid points

    u0 = exp.(sin.(2x))                             # Initial condition
    Fu0 = fft(u0)                                   # Frequency domain initial condition
     
    V = (1 .+ sin.(x).^2).^(-1)                     # Potential function

    L = DiagonalOperator(-im * k.^2)                # Linear part
    NL(u, p, t) = begin                             # Nonlinear part
        IFu = ifft(u)
        -im * fft((V + (abs.(IFu)).^2) .* IFu)
    end
    
    problem = SplitODEProblem(L, NL, Fu0, tspan)    # Define split u' = Lu + Nu problem
    sol = solve(problem, method; kwargs...)         # Solve the system

    u = ifft.(sol.u)                                # Convert solution back to real space

    return (sol.t, x, u)
end

function schrodinger_heatmap(N = 128, tspan = (0, 1), method = ETDRK4(autodiff = false); kwargs...)
    t, x, u = schrodinger(N, tspan, method; kwargs...)

    U = hcat(u...)

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax2 = Axis(fig[2,1],
              xlabel = "x",
              ylabel = "t",
              title = L"\mathrm{Re}\, \Psi"
             )
    ax3 = Axis(fig[2,2],
              xlabel = "x",
              ylabel = "t",
              title = L"\mathrm{Im}\, \Psi"
             )
    ax1 = Axis(fig[1,:],
              xlabel = "x",
              ylabel = "t",
              title = L"|\Psi|^2"
             )

    plt1 = heatmap!(ax1, x, t, abs.(U).^2, interpolate = true)
    plt2 = heatmap!(ax2, x, t, real(U), interpolate = true)
    plt3 = heatmap!(ax3, x, t, imag(U), interpolate = true)

    colsize!(fig.layout, 1, Aspect(1,1))
    colsize!(fig.layout, 2, Aspect(1,1))
    resize_to_layout!(fig)

    return fig
end

function schrodinger_animation(filename = "schrodinger", N = 128, tspan = (0,1), method = ETDRK4(autodiff = false); kwargs...)
    t, x, u = schrodinger(N, tspan, method; kwargs...)

    frame = Observable(1)
    U = @lift(u[$frame])
    Uabs  = @lift(abs.(u[$frame]).^2)
    Ureal = @lift(real(u[$frame]))
    Uimag = @lift(imag(u[$frame]))
    u_matrix = hcat(u...)

    abslims  = (minimum(abs.(u_matrix).^2), maximum(abs.(u_matrix)).^2)
    reallims = (minimum(real(u_matrix)), maximum(real(u_matrix)))
    imaglims = (minimum(imag(u_matrix)), maximum(imag(u_matrix)))

    fig = Figure()

    ax12 = Axis(fig[2,1],   xlabel = L"x", ylabel = "Time", title = L"\mathrm{Re}\,\psi", limits = (nothing, reallims))
    ax13 = Axis(fig[2,2],   xlabel = L"x", ylabel = "Time", title = L"\mathrm{Im}\,\psi", limits = (nothing, imaglims))
    ax11 = Axis(fig[1,1:2], xlabel = L"x", ylabel = "Time", limits = (nothing, abslims), 
                title = @lift("t = $(round(t[$frame], digits = 2))"))

    ax22 = Axis(fig[2,3],   xlabel = L"x", ylabel = "Time", title = L"\mathrm{Re}\,\psi")
    ax23 = Axis(fig[2,4],   xlabel = L"x", ylabel = "Time", title = L"\mathrm{Im}\,\psi")
    ax21 = Axis(fig[1,3:4], xlabel = L"x", ylabel = "Time", 
                title = @lift("t = $(round(t[$frame], digits = 2))"))

    lines!(ax11, x, Uabs,  color = :orange)
    lines!(ax12, x, Ureal, color = :orange)
    lines!(ax13, x, Uimag, color = :orange)

    heatmap!(ax21, x, t, abs.(u_matrix).^2, interpolate = true)
    heatmap!(ax22, x, t, real(u_matrix),    interpolate = true)
    heatmap!(ax23, x, t, imag(u_matrix),    interpolate = true)

    lines!(ax21, x, @lift(repeat([t[$frame]], length(x))), color = :orange, linewidth = 2)
    lines!(ax22, x, @lift(repeat([t[$frame]], length(x))), color = :orange, linewidth = 2)
    lines!(ax23, x, @lift(repeat([t[$frame]], length(x))), color = :orange, linewidth = 2)

    colsize!(fig.layout, 1, Aspect(1,1))
    colsize!(fig.layout, 2, Aspect(1,1))
    colsize!(fig.layout, 3, Aspect(1,1))
    colsize!(fig.layout, 4, Aspect(1,1))
    resize_to_layout!(fig)

    framerate = 60
    timestamps = 1:length(t)

    record(fig, filename*".mp4", timestamps; framerate = framerate) do i
        frame[] = i
    end
end
