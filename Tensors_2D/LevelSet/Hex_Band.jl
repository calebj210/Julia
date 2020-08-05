### Hex Band Generation Test

# Needed packages
using Plots
include("LevelSet.jl")

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Test function
function hexTest(N, n, r, c)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c;y(r,t)'.+c y(r,t)'.-c]

    # Compute hex band
    nodes = hexBand(zeroSet, N=10N, n=n)

    # Display number of nodes
    display(size(nodes,2))

    # Plot hex band
    plotA = scatter(nodes[1,:],nodes[2,:],
                    legend = false,
                    colorbar = true,
                    markeralpha = 0.75,
                    markersize = 3,
                    ratio = 1)
    display(plotA)
end

# Run test functions
hexTest(100, 10, 1, .85)
