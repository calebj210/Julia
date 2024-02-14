#=
# Test cases for region detection
#
# Author: Caleb Jacobs
# DLM: January 26, 2024
=#

include("RegionDetection.jl")
include("BackgroundNodes.jl")
using LinearAlgebra
using Plots

function cartCircTest(circN, N, M, n)
    Γ = circ(circN)                         # Circle nodes
    Ω = cartGrid(N, M, (-2, 2), (-2, 2))    # Background Cartesian nodes

    io = classifyNodes(Γ, Ω, n)             # Classification of nodes 

    trueIO = zeros(size(Ω, 2))              # True classification
    for i ∈ 1:size(Ω, 2)
        trueIO[i] = sign(1 - norm(Ω[:, i]))
    end

    p1 = scatter(Ω[1, :], Ω[2, :],
                 mz = io,
                 c = :viridis)
    p2 = scatter!(Γ[1, :], Γ[2, :],
                  ms = 5,
                  mα = 0.75,
                  legend = false,
                  aspect_ratio = :equal,
                  xlims = (-2.2, 2.2),
                  ylims = (-2.2, 2.2))

    display(p2)

    return [io trueIO (io - trueIO)]
end


function doubleCircTest(circN, N, M, n)
    Γ = circ(circN)                         # Circle nodes
    Ω = cartGrid(N, M, (-2, 2), (-2, 2))    # Background Cartesian nodes

    io = classifyNodes(Γ, Ω, n)             # Classification of nodes 

    trueIO = zeros(size(Ω, 2))              # True classification
    for i ∈ 1:size(Ω, 2)
        trueIO[i] = sign(1 - norm(Ω[:, i]))
    end

    p1 = scatter(Ω[1, :], Ω[2, :],
                 mz = io,
                 c = :viridis)
    p2 = scatter!(Γ[1, :], Γ[2, :],
                  ms = 5,
                  mα = 0.75,
                  legend = false,
                  aspect_ratio = :equal,
                  xlims = (-2.2, 2.2),
                  ylims = (-2.2, 2.2))

    display(p2)

    return [io trueIO (io - trueIO)]
end
