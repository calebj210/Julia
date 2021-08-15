using Plots
include("/home/merlin/Documents/Julia/Tensors_2D/Packages/Analysis.jl/src/Analysis.jl")
using ColorSchemes
using LaTeXStrings
using SparseArrays
include("LevelSet.jl")

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Radius over time
rt(r0,t) = sqrt(r0^2 + 2t)

##
# Generate band of helping nodes around the surface
#
function helpingNodes(zeroSet, norms, width)
    nodes = hcat(zeroSet + width * norms,
                 zeroSet - width * norms)

    return nodes
end

##
# Compute constrained interpolation weights
#
function getλ(nodes, norms, N, m)
    # "Collocation" matrix
    A  = [   ϕ(nodes[:, i], nodes[:, j], m)    for i ∈ 1:N, j ∈ 1:3N]
    Ax = [ϕ_xi(nodes[:, i], nodes[:, j], m, 1) for i ∈ 1:N, j ∈ 1:3N]
    Ay = [ϕ_xi(nodes[:, i], nodes[:, j], m, 2) for i ∈ 1:N, j ∈ 1:3N]

    # Constraint vector
    b  = zeros(N)
    bx = norms[1, :]
    by = norms[2, :]

    λ = [A; Ax; Ay] \ [b; bx; by]       # Constrained interpolation weights

    return λ
end

##
# Compute the next time step of the surface
#
# Updates zeroSet and norms
#
function nextStep!(zeroSet, norms, Δt, width, m)
    N = size(zeroSet, 2)                            # Number of surface nodes

    backNodes = helpingNodes(zeroSet, norms, width) # Helping background nodes

    nodes = hcat(zeroSet, backNodes)                # Combined node sets

    λ = getλ(nodes, norms, N, m)                    # Constrained interpolation weights

    # Gradient of distance function
    ∇f = hcat([S_xi(nodes, nodes[:, i], λ, 1, m = m) for i ∈ 1:3N],
              [S_xi(nodes, nodes[:, i], λ, 2, m = m) for i ∈ 1:3N])'

    # Distance function
    f = [(S(nodes, nodes[:, i], λ, m = m) - Δt * norm(∇f[:,i])) for i ∈ 1:3N]

    plotA = scatter(nodes[1, :], nodes[2,:],
                    marker_z    = f,
                    clims       = (-1,1),
                    c           = :viridis,
                    legend      = true,
                    colorbar    = true,
                    markeralpha = 0.75,
                    markersize  = 3,
                    ratio       = 1)
    display(plotA)
    sleep(5)

    # Run Coul-Newton to update surface
    zeroSet[:,:] = coulNewtonAdapt(zeroSet, nodes, f,
                                    n=11, m=m, o=2,
                                    maxIts=500, μ=1.25,
                                    Δt=10^(-5), ε=10^(-13))

    # Update surface normal vectors
    norms[:, :] = hcat([S_xi(nodes, zeroSet[:, i], λ, 1, m = m) for i ∈ 1:N],
                       [S_xi(nodes, zeroSet[:, i], λ, 2, m = m) for i ∈ 1:N])'
    for i ∈ 1:N
        norms[:, i] /= norm(norms[:, i])
    end

    return
end

function timeSolve(N; n, m, o, r, c, width, Δt, tf)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c; y(r,t)'.+c y(r,t)'.-c]

    N = 2N                  # Update number of nodes

    # Normal directions to the zero-set
    norms = zeros(2, N)
    for i ∈ 1:N
        tmp = [0,0]
        if zeroSet[2,i] >= -zeroSet[1,i]
            tmp = zeroSet[:, i] - [c,c]
        else
            tmp = zeroSet[:, i] + [c,c]
        end
        norms[:, i] = tmp / norm(tmp)
    end

    # Plot solution with errors overlayed
    for i = 0:10
        plotA = Analysis.vectorPlot(zeroSet, norms)
        display(plotA)
        sleep(5)

        nextStep!(zeroSet, norms, Δt, width, m)
    end
end

timeSolve(50, n=15, m=5, o=2, r=1, c=0.8, width=0.1, Δt = 10^(-3), tf = 0.05)
