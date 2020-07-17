### Reinitialization tester

# Needed packages
using Plots
using Analysis
using ColorSchemes
include("LevelSet.jl")

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Reinitialization test
function reinitTest(N; n, m, o, r)
    # Coordinate discretization
    zN = round(Int, N)
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

    # Reinitialize ambient nodes
    f = reinit(Nodes, oldNodes=nodes, F=F, zeroSet = zeroSet)
    

    # Compute new zero-set
    t = range(0,2*π-2*π/50, length = 50)
    zeroSet = [x(r,t)'; y(r,t)']

    zeroSet = coulNewtonAdapt(zeroSet, Nodes, f, n=n, m=m, o=o, maxIts=200, μ=2,
                              Δt=0.00001)

    # Compute magnitudes of new zero-set
    mags = zeros(50)
    for i∈1:50
        mags[i] = norm(zeroSet[:,i])
    end

    trues = fill(r, 50)

    tmp,errs = scalError(trues, mags)
    display(tmp)
    
    # # Construct distance function for ambient nodes
    # F = zeros(size(nodes,2))
    # for i∈eachindex(F)
    #     F[i] = norm(nodes[:,i]) - r
    # end
    
    # plotA = scatter(Nodes[1,:],Nodes[2,:],
    #                 marker_z = f,
    #                 c = :viridis,
    #                 legend = false,
    #                 colorbar = true,
    #                 colorbar_title = "Distance Function",
    #                 markeralpha = 0.75,
    #                 markerstrokewidth = 0,
    #                 ratio = 1,
    #                 markersize = 3,
    #                 xlims = (-1.25r,1.25r),
    #                 ylims = (-1.25r,1.25r),
    #                 clims = (-0.25,0.25))
    # t2 = range(0,2*π, length = 10N)
    # plot!(zeroSet[1,:],zeroSet[2,:],
    #       linewidth = 2,
    #       linealpha = 0.75,
    #       c = :black)
    plotA = errPlot(Nodes, scalError(ftrue, f)[2])
    # plotA = scatter(zeroSet[1,:],zeroSet[2,:],
    #                 marker_z = errs,
    #                 c = :viridis,
    #                 legend = false,
    #                 colorbar = true,
    #                 markeralpha = 0.75,
    #                 markerstrokewidth = 0,
    #                 ratio = 1,
    #                 markersize = 3,
    #                 xlims = (-1.25r,1.25r),
    #                 ylims = (-1.25r,1.25r))
    # plotA = errPlot(zeroSet, errs)
    display(plotA)
end


reinitTest(1000, n=7, m=5, o=1, r=1)
