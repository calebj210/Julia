### Level-Set Initialization Test

# Packages and functions to include
include("LevelSet.jl")
using Plots
using Analysis
using LaTeXStrings

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Test function for initialization on single circle
function oneInit(N; n, m, o, r, maxIts, ε, Δt)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)

    # Compute surface nodes
    zeroSet = [x(r,t)';y(r,t)']

    # Compute nearest neighbors among the surface
    kdtree = KDTree(zeroSet)
    idx = knn(kdtree, zeroSet, n, true)[1]

    # Compute approximate normals to surface
    norms = approxNormals(zeroSet, idx)
    norms = orVecs(zeroSet, norms, idx)

    # Re-orient surface if needed
    if norms[:,1]⋅[1,0] < 0
        norms *= -1
    end
    
    # Construct ambient nodes
    nodes = hexGen(10N, minx=-2r, maxx=2r, miny=-2r, maxy=2r)
    N = size(nodes,2)
    
    # Compute distance function for ambient nodes
    F = initialize(zeroSet, nodes, norms, n=n, m=m, o=o, maxIts=maxIts, ε=ε, Δt=Δt)

    # Compute true solution
    trues = zeros(N)
    for i∈1:N
        trues[i] = norm(nodes[:,i])-r
    end

    # Compute errors
    maxErr, err = scalError(trues, F)

    display(maxErr)
    
    # Plot functions and errors
    plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                    c = :steelblue,
                    label = "Surface Nodes",
                    ratio = 1,
                    dpi = 300,
                    markersize = 2,
                    xlims = (-2r,2r),
                    ylims = (-2r,2r))
    scatter!(nodes[1,:],nodes[2,:],
             markershape = :cross,
             label = "Background Nodes",
             c = :orange,
             markersize = 2)
                    
    plotB = scatter(nodes[1,:],nodes[2,:],
                    marker_z = F,
                    c = :viridis,
                    colorbar = true,
                    label = "Background Nodes",
                    colorbar_title = "Distance Function",
                    markeralpha = 0.75,
                    markerstrokewidth = 0,
                    ratio = 1,
                    markersize = 3,
                    dpi = 300,
                    xlims = (-2r,2r),
                    ylims = (-2r,2r))
    scatter!(zeroSet[1,:],zeroSet[2,:],
                    c = :steelblue,
                    label = "Surface Nodes",
                    ratio = 1,
                    markersize = 2)
    
    plotC = errPlot(nodes, err)
    plotC = plot(plotC,
                 colorbar_title = L"Log_{10}~Error",
                 legend = true,
                 label = "Background Nodes",
                 xlims = (-2r,2r),
                 ylims = (-2r,2r),
                 dpi = 300)
    scatter!(zeroSet[1,:],zeroSet[2,:],
             c = :steelblue,
             label = "Surface Nodes",
             ratio = 1,
             markersize = 2)
    
    display(plotA)
    readline()
    display(plotB)
    readline()
    display(plotC)
    readline()
    
    # Export plots
    savefig(plotA, "fig_1a.png")
    savefig(plotB, "fig_1b.png")
    savefig(plotC, "fig_1c.png")
end

# Initialization of two circles
function doubleInit(N; n, m, o, r, c, maxIts, ε, Δt)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c;y(r,t)'.+c y(r,t)'.-c]

    # Compute nearest neighbors among the surface
    kdtree = KDTree(zeroSet)
    idx = knn(kdtree, zeroSet, n, true)[1]

    # Compute approximate normals to surface
    norms = approxNormals(zeroSet, idx)
    norms = orVecs(zeroSet, norms, idx)

    # Re-orient surface if needed
    if norms[:,1]⋅[1,0] < 0
        norms[:,1:N] *= -1
    end
    if norms[:,N+1]⋅[1,0] < 0
        norms[:,N+1:end] *= -1
    end
    
    # Construct ambient nodes
    nodes = hexGen(10N, minx=-(c+2r), maxx=(c+2r), miny=-(c+2r), maxy=(c+2r))
    N = size(nodes,2)
    
    # Compute distance function for ambient nodes
    F = initialize(zeroSet, nodes, norms, n=n, m=m, o=o, maxIts=maxIts, ε=ε, Δt=Δt)

    # Compute true solution
    trues = zeros(N)
    for i∈1:N
        if nodes[2,i] >= -nodes[1,i]
            trues[i] = norm(nodes[:,i] - [c,c]) - r
        else
            trues[i] = norm(nodes[:,i] + [c,c]) - r
        end
    end

    # Compute errors
    maxErr, err = scalError(trues, F)

    display(maxErr)
    
    # Plot functions and errors
    plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                    c = :steelblue,
                    label = "Surface Nodes",
                    ratio = 1,
                    dpi = 300,
                    markersize = 2,
                    xlims = (-(c+2r),(c+2r)),
                    ylims = (-(c+2r),(c+2r)))
    scatter!(nodes[1,:],nodes[2,:],
             markershape = :cross,
             label = "Background Nodes",
             c = :orange,
             markersize = 2)
                    
    plotB = scatter(nodes[1,:],nodes[2,:],
                    marker_z = F,
                    c = :viridis,
                    colorbar = true,
                    label = "Background Nodes",
                    colorbar_title = "Distance Function",
                    markeralpha = 0.75,
                    markerstrokewidth = 0,
                    ratio = 1,
                    markersize = 3,
                    dpi = 300,
                    xlims = (-(c+2r),(c+2r)),
                    ylims = (-(c+2r),(c+2r)))
    scatter!(zeroSet[1,:],zeroSet[2,:],
                    c = :steelblue,
                    label = "Surface Nodes",
                    ratio = 1,
                    markersize = 2)
    
    plotC = errPlot(nodes, err)
    plotC = plot(plotC,
                 colorbar_title = L"Log_{10}~Error",
                 legend = true,
                 label = "Background Nodes",
                 xlims = (-(c+2r),(c+2r)),
                 ylims = (-(c+2r),(c+2r)),
                 dpi = 300)
    scatter!(zeroSet[1,:],zeroSet[2,:],
             c = :steelblue,
             label = "Surface Nodes",
             ratio = 1,
             markersize = 2)
    
    display(plotA)
    readline()
    display(plotB)
    readline()
    display(plotC)
    readline()

    # Export plots
    savefig(plotA, "fig_2a.png")
    savefig(plotB, "fig_2b.png")
    savefig(plotC, "fig_2c.png")
end

# Run test functions
oneInit(100, n=11, m=5, o=9, r=1, maxIts=500, ε=10^(-15), Δt=0.5)
doubleInit(100, n=11, m=5, o=9, r=1, c=1, maxIts=500, ε=10^(-15), Δt=0.3)
