#=
# Surface growth based on function density
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

function fibBarbell(N = 1000)
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
            size  = 40,
            color = u,
            colorscale = "Viridis",
#             cmin = -1,
#             cmax = 1,
            colorbar = attr(title="Concentration"),
            opacity = 1),
        ))
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
function densityGrowth(U, p, t)
    n,m,o = round.(Int, p[1:3])

    N = size(U, 2)

    u = U[1, :]
    nodes = U[2:4, :]
    
    c = getCommons(nodes, n, m, o)

    n⃗  = getNormals(nodes, c)
    Δu = computeΔ(nodes, u, c)
    H  = getH(nodes, c)
    divVn⃗ = computeSurfDiv(nodes, u, c, 1, 1)
#     divVn⃗ = computeSurfDivU(nodes, u, c)

    dU = zeros(4, N)
    dU[1, :] = 0.1Δu - divVn⃗
    for i ∈ 1:N
        dU[2:4, i] = (H[i] + u[i]) * n⃗[:, i]
    end
    
    return dU
end

function test1(T::Float64 = 15.0; N = 1000, n = 11, m = 3, o = 2)
    Γ = fibSphere(N)                        # Initial surface node-set

#     u = sum(abs, Γ', dims = 2) .- 1         # Initial u function
#     u = 0.9cos.(2π * Γ[3, :])               # Initial u function
#     u = 0.5cos.(2π * Γ[3, :]) - 0.5cos.(2π * Γ[1, :]) + 0.5cos.(2π * Γ[2, :])               # Initial u function
    u = 0.9cos.(π * Γ[3, :]) + 4 * Γ[3,:] .^ 10              # Initial u function
    
    U0 = [u'; Γ]
    tspan = (0.0, T)
    p = [n, m, o]
    
    prob = ODEProblem(densityGrowth, U0, tspan, p)
    sol  = solve(prob, 
                 alg_hints = [:stiff], 
                 reltol = 1e-2, 
                 abstol = 1e-2, 
                 save_everystep = false,
#                  saveat = 0.5,
                 progress = true,
                 progress_steps = 1)

    for i ∈ 1:length(sol)
        display(plotSurface(sol[i][2:4, :], sol[i][1, :]))
    end
        


    return sol
end

function test2(T::Float64 = 15.0; N = 1000, n = 11, m = 3, o = 2)
    Γ = fibBarbell(N)                        # Initial surface node-set

#     u = sum(abs, Γ', dims = 2) .- 1         # Initial u function
#     u = ℯ .^ (-(2Γ[3, :]).^2) .- 0.2                # Initial u function
#     u = 0.5cos.(2π * Γ[1, :])               # Initial u function
    u = 0.9cos.(π * Γ[3, :]) + 4 * Γ[3,:] .^ 10              # Initial u function
    
    U0 = [u'; Γ]
    tspan = (0.0, T)
    p = [n, m, o]
    
    prob = ODEProblem(densityGrowth, U0, tspan, p)
    sol  = solve(prob, 
                 alg_hints = [:stiff], 
                 reltol = 1e-2, 
                 abstol = 1e-2, 
                 save_everystep = false,
                 progress = true,
                 progress_steps = 1)

    display(plotSurface(sol[1][2:4, :], sol[1][1, :]))
    display(plotSurface(sol[end][2:4, :], sol[end][1, :]))

    return sol
end
