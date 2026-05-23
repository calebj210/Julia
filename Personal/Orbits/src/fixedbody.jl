using LinearAlgebra
import OrdinaryDiffEq as ODE
using GLMakie

function orb(; x0 = [1.0, 0.0], v0 = [0.0, 1.0], m = 1, T = (0.0, 10), kwargs...)
    function F(ddu, du, u, ω, t)
        ddu .= -ω * u / norm(u)^3
    end

    prob = ODE.SecondOrderODEProblem(F, v0, x0, T, m)
    sol = ODE.solve(prob, ODE.Tsit5(); kwargs...)

    fig = Figure()
    
    ax = Axis(fig[1,1], aspect = DataAspect())

    lines!(ax, sol, idxs = (3,4))

    resize_to_layout!(fig)

    return fig
end
