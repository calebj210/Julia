function practice()
    f(u, p, t) = 1.01u
    u0 = 1 / 2
    tspan = (0.0, 1.0)
    problem = ODEProblem(f, u0, tspan)
    sol = solve(problem, Tsit5(), reltol = 1e-16, abstol = 1e-16)

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax  = Axis(fig[1,1])
    ax.xlabel = "Time (t)"
    ax.ylabel = L"$u(t)$ (in $\mu$m)"
    ax.title  = "Solution to the linear ODE with a thick line"

    lines!(ax, sol, linewidth = 5, color = :blue, label = "My Thick Line!")
    lines!(sol.t, t -> .5exp(1.01t), linewidth = 3, linestyle = :dash, color = :orange, label = "True Solution")
    axislegend(ax, ax, position = :lt)

    ax2 = Axis(fig[2,1])
    lines!(ax2, sol.t, t -> abs((sol(t) - .5exp(1.01t)) / sol(t)), label = "Error")
    axislegend(ax2, ax2, position = :lt)

    return (fig, sol)
end

function lorenz()
    function lorenz!(du, u, p, t)
        du[1] = 10.0 * (u[2] - u[1])
        du[2] = u[1] * (28.0 - u[3]) - u[2]
        du[3] = u[1] * u[2] - (8 / 3) * u[3]
    end

    u0 = [1.0; 0.0; 0.0]
    tspan = (0.0, 100.0)
    prob = ODEProblem(lorenz!, u0, tspan)
    sol = solve(prob)

    fig = Figure()
    ax = Axis3(fig[1,1], aspect = (1,1,1))

    t = range(0,100,10000)
    lines!(ax, sol.(t, idxs = 1), sol.(t, idxs = 2), sol.(t, idxs = 3), color = t, colormap = :viridis)

    return (fig, sol)
end


