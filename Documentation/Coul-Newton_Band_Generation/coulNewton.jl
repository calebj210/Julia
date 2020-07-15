### Coulomb-Newton Method Tester

# Needed packages
using Plots
using ColorSchemes
include("LevelSet.jl")

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Narrow band test
function narrowBand(N; n, m, o, r)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Construct ambient nodes
    Nodes = hexGen(N, minx=-2r, maxx=2r, miny=-2r, maxy=2r)

    # Compute number of nodes in ambient node set
    NN = size(Nodes,2)
    
    # Construct distance function for ambient nodes
    F = zeros(NN)
    for i∈1:NN
        F[i] = norm(Nodes[:,i]) - r
    end

    # Construct narrow band
    nodes,F,NN = coulNewtonBand(Nodes, F, n=n, m=m, o=o, maxIts=500, μ=2,
                            η=0.01, Δt=0.00001, ε=0*10^(-3))

    # Construct distance function for ambient nodes
    F = zeros(size(nodes,2))
    for i∈eachindex(F)
        F[i] = norm(nodes[:,i]) - r
    end
    
    plotA = scatter(nodes[1,:],nodes[2,:],
                    marker_z = F,
                    c = :viridis,
                    legend = false,
                    colorbar = true,
                    colorbar_title = "Distance Function",
                    markeralpha = 0.75,
                    markerstrokewidth = 0,
                    ratio = 1,
                    markersize = 1.5,
                    xlims = (-1.25r,1.25r),
                    ylims = (-1.25r,1.25r))
    t2 = range(0,2*π, length = 10N)
    plot!(x(r,t2),y(r,t2),
          linewidth = 2,
          linealpha = 0.75,
          c = :black,
          dpi = 300)
    savefig(plotA, "fig_c.png")
    display(plotA)
end

# Narrow band test
function surfaceAdapt(N; n, m, o, r)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Create zero-level-set
    zeroSet = [x(r,t)';y(r,t)'] + 0.1r*(2*rand(Float64, (2,N)) .- 1)

    # Plot initial zero-set guess
    display(scatter(zeroSet[1,:],zeroSet[2,:],
                    ratio = 1))
    sleep(2)
    
    # Construct ambient nodes
    Nodes = hexGen(60N, minx=-2r, maxx=2r, miny=-2r, maxy=2r)

    # Compute number of nodes in ambient node set
    NN = size(Nodes,2)
    
    # Construct distance function for ambient nodes
    F = zeros(NN)
    for i∈1:NN
        F[i] = norm(Nodes[:,i]) - r
    end

    # Adapt zeroSet
    nodes = coulNewtonAdapt(zeroSet, Nodes, F, n=n, m=m, o=o, maxIts=100, μ=10, Δt=0.00001, ε=1*10^(-3))

    # Plot new zero set
    plotA = scatter(nodes[1,:],nodes[2,:],
                    legend = false,
                    ratio = 1,
                    markersize = 3)
    display(plotA)
end

# Growing circle test
# Narrow band test
function growingCircle(N; n, m, o, r, Δt, tf)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Create zero-level-set
    zeroSet = [x(r,t)';y(r,t)'] # + 0.2r*(2*rand(Float64, (2,N)) .- 1)
    
    # Construct ambient nodes
    nodes = hexGen(5N, minx=-1.5r, maxx=1.5r, miny=-1.5r, maxy=1.5r)

    # Compute number of nodes in ambient node set
    N = size(nodes,2)
    
    # Construct distance function for ambient nodes
    F = zeros(N)
    for i∈1:N
        F[i] = norm(nodes[:,i]) - r
    end

    # Construct narrow band
    nodes,F,N = coulNewton(zeroSet, nodes, F, n=n, m=m, o=o, maxIts=200, μ=3,
                           η=0.005, Δt=0.00001, ε=1*10^(-3))

    for i∈1:N
        F[i] = norm(nodes[:,i]) - r
    end

    # Compute total size of nodes
    NN = size(nodes,2)

    # Plot initial level-set function
    plotA = scatter(nodes[1,:],nodes[2,:],
                    marker_z = F,
                    legend = false,
                    colorbar = true,
                    markeralpha = 0.75,
                    markerstrokewidth = 0,
                    ratio = 1,
                    markersize = 1.5)
    display(plotA)
    sleep(2)

    # Discretize first order derivatives
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)
    
    # Begin evolution in time
    time = Δt:Δt:tf
    for i∈time
        readline()
        # Compute needed derivatives
        ϕx = Dx*F
        ϕy = Dy*F

        # Compute ∇ϕ
        ∇ϕ = [ϕx';ϕy']

        # Compute the normal direction: ∇ϕ/||∇ϕ||
        norm∇ϕ = zeros(2,NN)
        for j∈1:N
            norm∇ϕ[:,j] = ∇ϕ[:,j]/norm(∇ϕ[:,j])
        end

        # Compute curvature or the divergence of the normal direction
        κ = abs.(Dx*norm∇ϕ[1,:] + Dy*norm∇ϕ[2,:])

        # Compute (normal) ⋅ ∇ϕ
        nDot∇ϕ = zeros(NN)
        for j∈1:N
            nDot∇ϕ[j] = norm∇ϕ[:,j]⋅∇ϕ[:,j]
        end
        
        # Compute next time step
        F -= κ .* nDot∇ϕ * Δt

        # Find new zero-level-set
        # zeroSet = nodes[:,N+1:end]
        # zeroSet = newtonSolve(zeroSet, nodes, F, n=n, m=m, o=o,
        #                       ε=10^(-13), maxIts=5)
        
        # Plot zero level-set
        # plotA = scatter(nodes[1,:],nodes[2,:],
        #                 ratio = 1,
        #                 markersize = 1,
        #                 xlims = (-2r,2r),
        #                 ylims = (-2r,2r))
        # plotA = scatter(zeroSet[1,:],zeroSet[2,:],
        #                 ratio = 1,
        #                 markersize = 1,
        #                 xlims = (-2r,2r),
        #                 ylims = (-2r,2r))
        # display(plotA)

        # Plot initial level-set function
        plotA = scatter(nodes[1,:],nodes[2,:],
                 marker_z = F,
                 legend = false,
                 colorbar = true,
                 markeralpha = 0.75,
                 markerstrokewidth = 0,
                 ratio = 1,
                 markersize = 1.5)
        display(plotA)
    end

    display(plotA)
end


narrowBand(1000, n=7, m=5, o=1, r=1)
# surfaceAdapt(30, n=7, m=5, o=1, r=1)
# growingCircle(25, n=7, m=5, o=1, r=1, Δt=10^(-6), tf=0.0001)
