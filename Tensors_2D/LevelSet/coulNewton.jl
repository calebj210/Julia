### Coulomb-Newton Method Tester

# Needed packages
using Plots
include("LevelSet.jl")
gr()

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Narrow band test
function narrowBand(N; n, m, o, r)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Create zero-level-set
    zeroSet = [x(r,t)';y(r,t)'] + 0.2r*(2*rand(Float64, (2,N)) .- 1)

    display(scatter(zeroSet[1,:],zeroSet[2,:],
                    ratio = 1))
    sleep(2)
    
    # Construct ambient nodes
    nodes = hexGen(20N, minx=-2r, maxx=2r, miny=-2r, maxy=2r)

    # Compute number of nodes in ambient node set
    NN = size(nodes,2)
    
    # Construct distance function for ambient nodes
    F = zeros(NN)
    for i∈1:NN
        F[i] = norm(nodes[:,i]) - r
    end

    # Construct narrow band
    nodes,F,NN = coulNewton(zeroSet, nodes, F, n=n, m=m, o=o, maxIts=200, μ=0.75,
                       η=0.05, Δt=0.00001, ε=1*10^(-3))
    
    plotA = scatter(nodes[1,:],nodes[2,:],
                    marker_z = F,
                    c = :viridis,
                    legend = false,
                    colorbar = true,
                    markeralpha = 0.75,
                    markerstrokewidth = 0,
                    ratio = 1,
                    markersize = 1.5)
    # scatter!(nodes[1,NN+1:end],nodes[2,NN+1:end],
    #          markersize = 2)

    display(plotA)
end

narrowBand(25, n=7, m=5, o=1, r=.2)
