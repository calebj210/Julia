### Merging Circles Test

# Needed packages
using Plots
using Analysis
using ColorSchemes
using LaTeXStrings
include("LevelSet.jl")

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Test with hexagonal grid and no reinitialization
function hexNoReinit(N; n, m, o, r, c, maxIts=250, Δt, tf)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c;y(r,t)'.+c y(r,t)'.-c]

    # Generate node band
    nodes = hexGen(10N, minx=-(c+2r), maxx=(c+2r), miny=-(c+2r), maxy=(c+2r))
    N = size(nodes,2)

    # Initialize node band
    f = zeros(N)
    for i∈1:N
        if nodes[2,i] >= -nodes[1,i]
            f[i] = norm(nodes[:,i] - [c,c]) - r
        else
            f[i] = norm(nodes[:,i] + [c,c]) - r
        end
    end

    # Descretize required operators
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)

    # Initialize time variable
    time = Δt:Δt:tf

    # Initialize reinitialization counter
    reset = 0

    # Begin time evolution
    for i∈time
        # Compute needed derivatives
        ϕx = Dx*f
        ϕy = Dy*f

        # Compute ∇ϕ
        ∇ϕ = [ϕx';ϕy']

        # Compute the normal direction: ∇ϕ/||∇ϕ||
        norm∇ϕ = zeros(2,N)
        for j∈1:N
            norm∇ϕ[:,j] = ∇ϕ[:,j]/norm(∇ϕ[:,j])
        end

        # Compute curvature or the divergence of the normal direction
        κ = abs.(Dx*norm∇ϕ[1,:] + Dy*norm∇ϕ[2,:])

        # Compute normal ⋅ ∇ϕ
        nDot∇ϕ = zeros(N)
        for j∈1:N
            nDot∇ϕ[j] = norm∇ϕ[:,j]⋅∇ϕ[:,j]
        end
        
        # Compute next time step
        f -= κ .* nDot∇ϕ * Δt

        # Check for solution plotting
        if reset >= 100
            # zeroSet = newtonSolve(zeroSet, nodes, f, n=n, m=m, o=o,
            #                       ε=10^(-10), maxIts=5)
            zeroSet = coulNewtonAdapt(zeroSet, nodes, f, n=n, m=m, o=o,
                                      maxIts=50, μ=2, Δt=10^(-5), ε=10^(-10))
            
            # Plot level-set function
            # plotA = scatter(nodes[1,:],nodes[2,:],
            #                 marker_z = f,
            #                 c = :viridis,
            #                 legend = false,
            #                 colorbar = true,
            #                 markeralpha = 0.75,
            #                 markersize = 3,
            #                 ratio = 1)
            
            plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                            ratio = 1,
                            legend = false,
                            xlims = (-(c+2r),c+2r),
                            ylims = (-(c+2r),c+2r),
                            markersize = 2)
            display(plotA)
            
            reset = 0
        end

        reset += 1
    end
end

# Test function on hexagonal grid with reinitialization
function hexWithReinit(N; n, m, o, r, c, maxIts=250, Δt, tf)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c;y(r,t)'.+c y(r,t)'.-c]

    # Generate node band
    nodes = hexGen(10N, minx=-(c+2r), maxx=(c+2r), miny=-(c+2r), maxy=(c+2r))
    N = size(nodes,2)

    # Initialize node band
    f = zeros(N)
    for i∈1:N
        if nodes[2,i] >= -nodes[1,i]
            f[i] = norm(nodes[:,i] - [c,c]) - r
        else
            f[i] = norm(nodes[:,i] + [c,c]) - r
        end
    end

    # Descretize required operators
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)

    # Initialize time variable
    time = Δt:Δt:tf

    # Initialize reinitialization counter
    reset = 0

    # Begin time evolution
    for i∈time
        # Check for reinitialization
        if reset % 100 == 0
            f = reinit(nodes, oldNodes=nodes, F=f, zeroSet=zeroSet, smooth=true, n=n, m=m, o=o)
        end
        
        # Compute needed derivatives
        ϕx = Dx*f
        ϕy = Dy*f

        # Compute ∇ϕ
        ∇ϕ = [ϕx';ϕy']

        # Compute the normal direction: ∇ϕ/||∇ϕ||
        norm∇ϕ = zeros(2,N)
        for j∈1:N
            norm∇ϕ[:,j] = ∇ϕ[:,j]/norm(∇ϕ[:,j])
        end

        # Compute curvature or the divergence of the normal direction
        κ = abs.(Dx*norm∇ϕ[1,:] + Dy*norm∇ϕ[2,:])

        # Compute normal ⋅ ∇ϕ
        nDot∇ϕ = zeros(N)
        for j∈1:N
            nDot∇ϕ[j] = norm∇ϕ[:,j]⋅∇ϕ[:,j]
        end
        
        # Compute next time step
        f -= κ .* nDot∇ϕ * Δt

        # Check for solution plotting
        if reset >= 100
            # zeroSet = newtonSolve(zeroSet, nodes, f, n=n, m=m, o=o,
            #                       ε=10^(-10), maxIts=5)
            zeroSet = coulNewtonAdapt(zeroSet, nodes, f, n=n, m=m, o=o, maxIts=50, μ=2, Δt=10^(-5), ε=10^(-10))
            
            # Plot level-set function
            # plotA = scatter(nodes[1,:],nodes[2,:],
            #                 marker_z = f,
            #                 c = :viridis,
            #                 legend = false,
            #                 colorbar = true,
            #                 markeralpha = 0.75,
            #                 markersize = 3,
            #                 ratio = 1)
            
            plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                            ratio = 1,
                            legend = false,
                            xlims = (-(c+2r),c+2r),
                            ylims = (-(c+2r),c+2r),
                            markersize = 2)
            sleep(0.1)
            display(plotA)

            reset = 0
        end

        reset += 1
    end
end

function narrowBandNoReinit(N; n, m, o, r, c, maxIts=250, Δt, tf)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c;y(r,t)'.+c y(r,t)'.-c]

    # Generate background nodes
    nodes = hexGen(10N, minx=-(c+2r), maxx=(c+2r), miny=-(c+2r), maxy=(c+2r))
    N = size(nodes,2)

    # Initialize node band
    f = zeros(N)
    for i∈1:N
        if nodes[2,i] >= -nodes[1,i]
            f[i] = norm(nodes[:,i] - [c,c]) - r
        else
            f[i] = norm(nodes[:,i] + [c,c]) - r
        end
    end
    
    # Generate node band
    nodes,f,N = coulNewtonBand(nodes, f, n=n, m=m, o=o, maxIts=500, μ=5, η=0.01, Δt=10^(-5), ε=0*10^(-3))

    plotA = scatter(nodes[1,:],nodes[2,:],
                    ratio = 1,
                    legend = false,
                    markersize=2)
    display(plotA)

    readline()
    
    # Descretize required operators
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)

    # Initialize time variable
    time = Δt:Δt:tf

    # Initialize reinitialization counter
    reset = 0

    # Begin time evolution
    for i∈time
        # Compute needed derivatives
        ϕx = Dx*f
        ϕy = Dy*f

        # Compute ∇ϕ
        ∇ϕ = [ϕx';ϕy']

        # Compute the normal direction: ∇ϕ/||∇ϕ||
        norm∇ϕ = zeros(2,N)
        for j∈1:N
            norm∇ϕ[:,j] = ∇ϕ[:,j]/norm(∇ϕ[:,j])
        end

        # Compute curvature or the divergence of the normal direction
        κ = abs.(Dx*norm∇ϕ[1,:] + Dy*norm∇ϕ[2,:])

        # Compute normal ⋅ ∇ϕ
        nDot∇ϕ = zeros(N)
        for j∈1:N
            nDot∇ϕ[j] = norm∇ϕ[:,j]⋅∇ϕ[:,j]
        end
        
        # Compute next time step
        f -= κ .* nDot∇ϕ * Δt

        # Check for solution plotting
        if reset >= 100
            zeroSet = newtonSolve(zeroSet, nodes, f, n=n, m=m, o=o,
                                  ε=10^(-10), maxIts=5)
            
            # Plot level-set function
            # plotA = scatter(nodes[1,:],nodes[2,:],
            #                 marker_z = f,
            #                 c = :viridis,
            #                 legend = false,
            #                 colorbar = true,
            #                 markeralpha = 0.75,
            #                 markersize = 3,
            #                 ratio = 1)
            
            plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                            ratio = 1,
                            legend = false,
                            xlims = (-(c+2r),c+2r),
                            ylims = (-(c+2r),c+2r),
                            markersize = 3)
            sleep(0.1)
            display(plotA)

            reset = 0
        end

        reset += 1
    end
end

# Test with hexagonal grid and no reinitialization
function hexBandNoReinit(N; n, m, o, r, c, maxIts=250, Δt, tf)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c;y(r,t)'.+c y(r,t)'.-c]

    # Generate node band
    nodes = hexBand(zeroSet, N=10N, n=10)
    N = size(nodes,2)

    # Plot hex band
    plotA = scatter(nodes[1,:],nodes[2,:],
                    legend = false,
                    colorbar = true,
                    markeralpha = 0.75,
                    markersize = 3,
                    ratio = 1)
    display(plotA)
    readline()

    # Initialize node band
    f = zeros(N)
    for i∈1:N
        if nodes[2,i] >= -nodes[1,i]
            f[i] = norm(nodes[:,i] - [c,c]) - r
        else
            f[i] = norm(nodes[:,i] + [c,c]) - r
        end
    end

    # Descretize required operators
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)

    # Initialize time variable
    time = Δt:Δt:tf

    # Initialize reinitialization counter
    reset = 0

    # Begin time evolution
    for i∈time
        # Compute needed derivatives
        ϕx = Dx*f
        ϕy = Dy*f

        # Compute ∇ϕ
        ∇ϕ = [ϕx';ϕy']

        # Compute the normal direction: ∇ϕ/||∇ϕ||
        norm∇ϕ = zeros(2,N)
        for j∈1:N
            norm∇ϕ[:,j] = ∇ϕ[:,j]/norm(∇ϕ[:,j])
        end

        # Compute curvature or the divergence of the normal direction
        κ = abs.(Dx*norm∇ϕ[1,:] + Dy*norm∇ϕ[2,:])

        # Compute normal ⋅ ∇ϕ
        nDot∇ϕ = zeros(N)
        for j∈1:N
            nDot∇ϕ[j] = norm∇ϕ[:,j]⋅∇ϕ[:,j]
        end
        
        # Compute next time step
        f -= κ .* nDot∇ϕ * Δt

        # Check for solution plotting
        if reset >= 200
            # zeroSet = newtonSolve(zeroSet, nodes, f, n=n, m=m, o=o,
            #                       ε=10^(-10), maxIts=5)
            zeroSet = coulNewtonAdapt(zeroSet, nodes, f, n=n, m=m, o=o, maxIts=50, μ=2, Δt=10^(-5), ε=10^(-10))
            
            # Plot level-set function
            # plotA = scatter(nodes[1,:],nodes[2,:],
            #                 marker_z = f,
            #                 c = :viridis,
            #                 legend = false,
            #                 colorbar = true,
            #                 markeralpha = 0.75,
            #                 markersize = 3,
            #                 ratio = 1)
            
            plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                            ratio = 1,
                            legend = false,
                            xlims = (-(c+2r),c+2r),
                            ylims = (-(c+2r),c+2r),
                            markersize = 2)
            sleep(0.1)
            display(plotA)

            reset = 0
        end

        reset += 1
    end
end

# Run test function
hexNoReinit(100, n=11, m=5, o=2, r=1, c=0.8, Δt=10^(-5), tf=0.1)

# hexWithReinit(100, n=11, m=5, o=2, r=1, c=0.8, Δt=10^(-4), tf=0.1)

# narrowBandNoReinit(100, n=11, m=5, o=2, r=1, c=0.8, Δt=10^(-4), tf=0.1)

# hexBandNoReinit(100, n=11, m=5, o=2, r=1, c=0.8, Δt=10^(-5), tf=0.1)
