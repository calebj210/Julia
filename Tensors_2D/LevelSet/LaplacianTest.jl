### Test for solving the heat equation in time

# Packages to include and use
include("LevelSet.jl")
using Plots

function lapTest(N,M;n,m,o,Δt,tf)
    # Form our grid
    x = -1:2/(N-1):1
    y = -1:2/(M-1):1
    nodes = [repeat(x,inner = (M,1))';
             repeat(y,outer = (N,1))']

    # Construct initial conditions
    f = zeros(N*M) .+ 5
    f[73] = 10
    # f[1:M] .= 0
    # f[M:M:N*M] .= 0
    # f[(N-1)*M+1:N*M] .= 0
    # f[1:M:N*M-1] .= 0

    # Compute discretized laplacian
    D = I-Δt*discretize∇²(nodes,n,m,o)

    # Plot initial temperature
    plotA = scatter(nodes[1,:],nodes[2,:],
                    marker_z = f,
                    colorbar = :right,
                    aspectratio = :equal,
                    legend = false,
                    markersize = 5,
                    clims = (0,10))
    display(plotA)
    sleep(1)
    
    # Evolve our PDE
    for t = Δt:Δt:tf
        # Solve for next step
        f = D\f

        # f[1:M] .= 0
        # f[M:M:N*M] .= 0
        # f[(N-1)*M+1:N*M] .= 0
        # f[1:M:N*M-1] .= 0

        # Plot temperature
        plotA = scatter(nodes[1,:],nodes[2,:],
                        marker_z = f,
                        colorbar = :right,
                        aspectratio = :equal,
                        legend = false,
                        markersize = 5,
                        clims = (0,10))
        display(plotA)
        sleep(0.01)
    end
end

lapTest(10,10, n=5, m=11, o=1,Δt = 10^(-2), tf = 1)
