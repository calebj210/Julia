### Testing the addition of nodes to our surface
using Plots
using Analysis
include("LevelSet.jl")

# Unit circle definition
# Unit circle
x(t,r) = r*cos.(t)
y(t,r) = r*sin.(t)


# Node band test
function nodesTest(N = 20, n=5, m=3, o=0; r=0.2)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Creating our node distribution
    nodes = [x(t,r)';y(t,r)']

    # Add node bands to node set and orient the surface
    (nodes,f) = generateNodeBand(nodes)
    if nodes[1,1] - nodes[1,N+1] < 0
        f *= -1
    end
    
    # Plot our new node set
    plotA = scatter(nodes[1,:],nodes[2,:],
                    marker_z = f,
                    colorbar = :right,
                    aspectratio = :equal,
                    legend = false,
                    markersize = 2)
    display(plotA)
end

# Normal test
function normsTest(N=20, n=5, m=3, o=0; r=0.2)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Creating our node distribution
    nodes = [x(t,r)';y(t,r)']

    # Add node bands to node set and orient the surface
    (nodes,f) = generateNodeBand(nodes)
    if nodes[1,1] - nodes[1,N+1] < 0
        f *= -1
    end

    # Find normals
    norms = findNormals(nodes, f, n, m, o)

    # Plot our normals
    plotA = vectorPlot(nodes[:,1:N],norms[:,1:N])
    display(plotA)

    # Errors
    err,normErrs = vecError(nodes[:,1:N],norms[:,1:N])

    display(err)
end

# Curvature test
function κTest(N=20, n=5, m=3, o=0; r=0.2)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Creating our node distribution
    nodes = [x(t,r)';y(t,r)']

    # Add node bands to node set and orient the surface
    (nodes,f) = generateNodeBand(nodes)
    if nodes[1,1] - nodes[1,N+1] < 0
        f *= -1
    end

    # Compute curvature
    κ = computeκ(nodes, f, n, m, o)

    # True solution
    trueκ = zeros(N)
    for i ∈ 1:N
        trueκ[i] = 1/norm(nodes[:,i])
    end
    
    # Compute the errors
    err,normErrs = scalError(trueκ, κ)
    
    # Display the error
    display(err)
    # display(κ)
end

# Growing circle test
function growTest(N=20, n=5, m=3, o=0; r=0.2, tf = 0.2, Δt = 10^(-4))
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Creating our node distribution
    nodes = [x(t,r)';y(t,r)']

    # Initialize our time variable
    time = Δt:Δt:tf
    
    # Begin surface growth
    for i∈time
        # Add node bands to node set and orient the surface
        (nodes,f) = generateNodeBand(nodes[:,1:N])
        if nodes[1,1] - nodes[1,N+1] < 0
            f *= -1
        end
        
        # Find normals
        norms = findNormals(nodes, f, n, m, o)
        
        # Compute curvature
        κ = computeκ(nodes, f, n, m, o)
        
        # Compute κn̄
        norms .*= [κ';κ']

        # Advance to next time step
        nodes += Δt*norms

        # Plot our growth
        plotA = scatter(nodes[1,1:N],nodes[2,1:N],
                        marker_z = f,
                        colorbar = :right,
                        aspectratio = :equal,
                        legend = false,
                        markersize = 2,
                        xlims = (-1,1),
                        ylims = (-1,1))
        # plotB = vectorPlot(nodes[:,1:N], norms[:,1:N])
        display(plotA)
        sleep(0.01)
    end
end

# Test for Newton's method
function newtonTest(N = 20; n=5, m=3, o=0, r=0.2)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Creating our node distribution
    nodes = [x(t,r)';y(t,r)']

    # Add node bands to node set and orient the surface
    (nodes,f) = generateNodeBand(nodes)
    if nodes[1,1] - nodes[1,N+1] < 0
        f *= -1
    end

    # Compute initial guess
    X = nodes[:,N+1:2N]

    # Compute zero level set
    X = newtonSolve(X, nodes, f, n=n, m=m, o=o, ε=10^(-13), maxIts = 100)

    # Compute magnitude of each node
    # for i∈1:N
    #     display(norm(X[:,i]))
    # end
    
    # # Plot our new node set
    plotA = scatter(nodes[1,N+1:end],nodes[2,N+1:end],
                    marker_z = f,
                    colorbar = :right,
                    aspectratio = :equal,
                    legend = false,
                    markersize = 1.5)
    scatter!(X[1,:],X[2,:],
             markersize = 1)
    display(plotA)
end

# nodesTest(25, 5, 3)
# normsTest(60, 5, 7, -1, r = 1)
# κTest(100, 11, 7, 2, r=0.2)
# growTest(60, 11, 7, 2, r=0.2, tf=0.01, Δt = 10^(-4))
newtonTest(25, n=5, m=5, o=-1, r=1)
