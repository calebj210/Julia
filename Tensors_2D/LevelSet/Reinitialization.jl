### Reinitialization tester

# Needed packages
using Plots
using Analysis
using ColorSchemes
using LaTeXStrings
include("LevelSet.jl")

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Reinitialization test
function reinitTest(N; n, m, o, r)
    # Coordinate discretization
    zN = round(Int, 10N)
    t = range(0,2*π-2*π/zN, length = zN)
    
    # Construct ambient nodes
    Nodes = hexGen(N, minx=-2r, maxx=2r, miny=-2r, maxy=2r)

    # Construct zero-set
    zeroSet = [x(r,t)'; y(r,t)']

    # Compute number of nodes in ambient node set
    N = size(Nodes,2)
    
    # Construct distance function for ambient nodes
    F = zeros(N)
    for i∈1:N
        F[i] = norm(Nodes[:,i]) - r
    end
    
    # Compute true f
    ftrue = zeros(N)
    for i∈1:N
        ftrue[i] = norm(Nodes[:,i]) - r
    end

    # Construct narrow band
    nodes,F,NN = coulNewtonBand(Nodes, F, n=n, m=m, o=o, maxIts=500, μ=2,
                                η=0.01, Δt=0.00001, ε=0*10^(-3))

    # Initialize narrow band
    F = zeros(size(nodes,2))
    for i∈eachindex(F)
        F[i] = norm(nodes[:,i]) - r
    end

    # Reinitialize ambient nodes
    f = reinit(Nodes, oldNodes=nodes, F=F, zeroSet = zeroSet, replace = false,
    n=n, m=m, o=o)
    
    # Compute new zero-set
    t = range(0,2*π-2*π/150, length = 150)
    zeroSet = [x(r,t)'; y(r,t)']

    zeroSet = coulNewtonAdapt(zeroSet, Nodes, f, n=n, m=m, o=o, maxIts=200, μ=2,
                              Δt=0.00001, ε=10^(-14))

    # Compute magnitudes of new zero-set
    mags = zeros(150)
    for i∈1:150
        mags[i] = norm(zeroSet[:,i])
    end

    trues = fill(r,150)

    tmp,errs = scalError(trues, mags)
    display(tmp)

    plotA = scatter(nodes[1,:],nodes[2,:],
                    marker_z = F,
                    c = :viridis,
                    legend = false,
                    colorbar = true,
                    title = "Narrow Band Distance Function",
                    colorbar_title = "Distance Function",
                    markeralpha = 0.75,
                    markerstrokewidth = 0,
                    ratio = 1,
                    markersize = 2.5,
                    dpi = 300)
    display(plotA)
    readline()
    # plotA = scatter(Nodes[1,:],Nodes[2,:],
    #                 marker_z = errs,
    #                 c = :viridis,
    #                 legend = false,
    #                 colorbar = true,
    #                 colorbar_title = "Distance Function",
    #                 markeralpha = 0.75,
    #                 markerstrokewidth = 0,
    #                 ratio = 1,
    #                 markersize = 3,
    #                 clims = (-2r,2r))
    # t2 = range(0,2*π, length = 10N)
    # plot!(zeroSet[1,:],zeroSet[2,:],
    #       linewidth = 2,
    #       linealpha = 0.75,
    #       c = :black)
    plotB = errPlot(Nodes, scalError(ftrue, f)[2])
    plotB = plot(plotB,
                 title = "Reinitialization Error",
                 colorbar_title = L"Log_{10}~Error",
                 dpi = 300)
    display(plotB)
    readline()
    plotC = scatter(zeroSet[1,:],zeroSet[2,:],
                    marker_z = log10.(errs),
                    c = :viridis,
                    legend = false,
                    colorbar = true,
                    colorbar_title = "Surface Error",
                    markeralpha = 0.75,
                    markerstrokewidth = 0,
                    ratio = 1,
                    markersize = 3,
                    xlims = (-1.25r,1.25r),
                    ylims = (-1.25r,1.25r))
    # plotC = errPlot(zeroSet, errs)
    plotC = plot(plotC,
                 title = "Surface Error",
                 colorbar_title = L"Log_{10}~Error",
                 dpi = 300)
    display(plotC)
    readline()

    savefig(plotA, "NarrowBand.png")
    savefig(plotB, "ReinitializedNodes.png")
    savefig(plotC, "Surface.png")
end



reinitTest(1000, n=7, m=5, o=1, r=1)
