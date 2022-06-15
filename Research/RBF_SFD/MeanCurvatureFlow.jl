#=
# Mean curvature flow simulation
# 
# Author: Caleb Jacobs
# DLM: 2022-06-14
=#

using DifferentialEquations
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

function fibBarbbell(N = 1000)
    ϕ = π * (3 - sqrt(5))
    
    Γ = zeros(3, N)
    for i ∈ 0:N - 1
        z = 1 - 2(i / (N - 1))
        r = (2z^2 + 0.2) * sqrt(1 - z^2)

        θ = ϕ * i

        x = r * cos(θ)
        y = r * sin(θ)

        Γ[:, i + 1] = [x, y, z]
    end

    return Γ
end

function plotSurface(Γ, u)
    p = plot(scatter(
        x    = Γ[1, :],
        y    = Γ[2, :],
        z    = Γ[3, :],
        type = "scatter3d",
        mode = "markers",
        marker = attr(
            size  = 10,
            color = u,
            colorscale = "Viridis",
            colorbar = attr(title="Concentration"),
            opacity = 1),
        ))
    return p
end

#= 
# RHS of tumor growth PDE
# 
# p = [n, m, o]
#   = [1, 2, 3]
=#
function meanCurvature(Γ, p, t)
    n, m, o = round.(Int, p[1:3])

    N = size(Γ, 2)
    
    c = getCommons(Γ, n, m, o)

    n⃗ = getNormals(nodes, c)
    H = getH(Γ, c)

    dΓ = zeros(3, N)
    for i ∈ 1:N
        dΓ[:, i] = H[i] * n⃗[:, i]
    end
    
    return dΓ
end

function test1(T::Float64 = 15.0; R = 1.0, N = 1000, n = 11, m = 3, o = 2)
    Γ0 = R * fibSphere(N)           # Initial surface node-set

    tspan = (0.0, T)
    p = [n, m, o]
    
    prob = ODEProblem(meanCurvature, Γ0, tspan, p)
    sol  = solve(prob, 
                 alg_hints = [:stiff], 
                 reltol = 1e-2, 
                 abstol = 1e-2, 
                 progress = true,
                 progress_steps = 1)

    display(plotSurface(sol[1], sol[1][3, :]))
    display(plotSurface(sol[end], sol[end][3, :]))

    return sol
end

function test2(T::Float64 = 15.0; N = 1000, n = 11, m = 3, o = 2)
    Γ0 = fibBarbbell(N)           # Initial surface node-set

    tspan = (0.0, T)
    p = [n, m, o]
    
    prob = ODEProblem(meanCurvature, Γ0, tspan, p)
    sol  = solve(prob, 
                 alg_hints = [:stiff],
                 reltol = 1e-2, 
                 abstol = 1e-2, 
                 progress = true,
                 progress_steps = 1)

    display(plotSurface(sol[1], sol[1][3, :]))
    display(plotSurface(sol[end], sol[end][3, :]))

    return sol
end
