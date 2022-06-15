#=
# Tumor growth simulation
# 
# Author: Caleb Jacobs
# DLM: 2022-06-14
=#

using DifferentialEquations
using Random
using PlotlyJS
using Printf
include("RBF_SFD.jl")

# Fibonacci sphere node generator
function fibSphere(N = 1000)
    ϕ = π * (3 - sqrt(5))       # Golden angle

    nodes = zeros(3, N)
    for i ∈ 0:N - 1
        y = 1 - 2(i / (N - 1))  # y-coordinate
        r = sqrt(1 - y^2)       # Radius

        θ = ϕ * i               # Angle

        x = r * cos(θ)          # x-coordinate
        z = r * sin(θ)          # y-coordinate

        nodes[:, i + 1] = [x, y, z]
    end

    return nodes
end
function fibTest(N = 1000)
    nodes = fibSphere(N)

    plot(scatter(
        x    = nodes[1, :],
        y    = nodes[2, :],
        z    = nodes[3, :],
        type = "scatter3d",
        mode = "markers",
        marker = attr(
            size = 2,
            colorscale = "Viridis",
            opacity = 0.9)))
end

function plotSurface(Γ, u; ttl = "")
    layout = Layout(
        title = attr(
            text = ttl,
            xanchor = "center",
            yanchor = "top"
        ),

        scene = attr(
            xaxis = attr(
                range = [-2, 2]
            ),

            yaxis = attr(
                range = [-2, 2]
            ),

            zaxis = attr(
                range = [-2, 2]
            )
        )

    )

    sctr = scatter(
        x    = Γ[1, :],
        y    = Γ[2, :],
        z    = Γ[3, :],
        type = "scatter3d",
        mode = "markers",
        marker = attr(
            size  = 20,
            color = u,
            colorscale = "Viridis",
            colorbar = attr(title=""),
            opacity = 1
        )
    )

    p = plot(sctr, layout)

    return p
end

#= 
# RHS of tumor growth PDE
# 
# U = [u; w; Γ]
#
# p = [n, m, o, a, b, d, γ, ε, δ]
#   = [1, 2, 3, 4, 5, 6, 7, 8, 9]
=#
function tumorGrowth(U, p, t)
    n,m,o = round.(Int, p[1:3])
    a,b,d,γ,ε,δ = p[4:end]

    N = size(U, 2)

    u = U[1, :]
    w = U[2, :]
    nodes = U[3:5, :]
    
    c = getCommons(nodes, n, m, o)

    n⃗  = getNormals(nodes, c)
    V  = computeNormVel(nodes, u, c, ε, δ)
    divVn⃗ = computeSurfDiv(nodes, u, c, ε, δ)
    Δu = computeΔ(nodes, u, c)
    Δw = computeΔ(nodes, w, c)
    F1 = [γ * (a - u[i] + u[i]^2 * w[i]) for i ∈ 1:N]
    F2 = [γ * (b -        u[i]^2 * w[i]) for i ∈ 1:N]

    dU = zeros(5, N)
    dU[1, :] =      Δu + F1 - u .* divVn⃗
    dU[2, :] =  d * Δw + F2 - w .* divVn⃗
    for i ∈ 1:N
        dU[3:5, i] = V[i] * n⃗[:, i]
    end
    
    return dU
end

function test1(T::Float64 = 1.0; ε = 0.01, δ = 0.4, N = 1000, n = 11, m = 3, o = 2)
#     a = 0.2
    a = 0.1
#     b = 1
    b = 0.9
    k = 6
#     d = 15
    d = 10
#     dc = 17.0056
#     γ = 2dc * k / (dc * (b - a) / (a + b) - (a + b)^2)
    γ = 100

    Γ = fibSphere(N)                        # Initial surface node-set

    Random.seed!(210)
#     u = repeat([a + b], N)  + 0.1 * Γ[1, :] # Initial u function
#     w = repeat([b / (a + b)^2], N)          # Initial w function
    u = 1 .+ 2Γ[1, :] .* Γ[2, :] .* Γ[3, :]
    w = 1 .+ 2Γ[1, :] .* Γ[2, :] .* Γ[3, :]

    U0 = [u'; w'; Γ]
    tspan = (0.0, T)
    p = [n, m, o, a, b, d, γ, ε, δ]
    
    prob = ODEProblem(tumorGrowth, U0, tspan, p)
    sol  = solve(prob, 
                 alg_hints = [:stiff], 
                 reltol = 1e-2, 
                 abstol = 1e-2, 
                 saveat = 0.5,
#                  save_everystep = false,
                 progress = true,
                 progress_steps = 1)

    for i ∈ 1:length(sol)
        display(plotSurface(sol[i][3:5, :], sol[i][1, :], ttl = "t = $(sol.t[i])"))
    end

    return sol
end


function test2(T::Float64 = 15.0; ε = 0.0, δ = 0.0, N = 1000, n = 11, m = 3, o = 2)
    a = 0.2
    b = 1
    k = 20
    d = 18
    dc = 17.0056
    γ = 2dc * k / (dc * (b - a) / (a + b) - (a + b)^2)

    Γ = fibSphere(N)                        # Initial surface node-set

    Random.seed!(210)
    u = repeat([a + b], N)  + 0.1 * Γ[1, :] # Initial u function
    w = repeat([b / (a + b)^2], N)          # Initial w function

    U0 = [u'; w'; Γ]
    tspan = (0.0, T)
    p = [n, m, o, a, b, d, γ, ε, δ]
    
    prob = ODEProblem(tumorGrowth, U0, tspan, p)
    sol  = solve(prob, 
                 alg_hints = [:stiff], 
                 reltol = 1e-2, 
                 abstol = 1e-2, 
                 save_everystep = false,
                 progress = true,
                 progress_steps = 1)

    display(plotSurface(sol[end][3:5, :], sol[end][1, :]))

    return sol
end
