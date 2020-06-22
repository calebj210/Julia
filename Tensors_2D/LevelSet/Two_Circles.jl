## Level narrow band level-set method for solving merging circles

# Packages to include and use
include("LevelSet.jl")
using Plots

# Suface functions
x(t) = 0.2*cos.(t)
y(t) = 0.2*sin.(t)

function grow(N; n, m, o, Δt, tf, initCnt)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Creating zero level-set
    zeroSet = zeros(2,2N)
    zeroSet[:,1:N] = [x(t)';y(t)'] .+ 0.2
    zeroSet[:,N+1:2N] = [x(t)';y(t)'] .- 0.2

    # Construct surface sets
    circ1 = [1:N...,2N+1:3N...,4N+1:5N...,6N+1:7N...,8N+1:9N...]
    circ2 = [N+1:2N...,3N+1:4N...,5N+1:6N...,7N+1:8N...,9N+1:10N...]
    
    # Plot full node set
    plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                     ratio = 1,
                     markersize = 2,
                    xlims = (-0.6,0.6),
                    ylims = (-0.6,0.6))
    display(plotA)
    sleep(0.25)
    
    # Initialize time variable
    time = Δt:Δt:tf

    # Initialize node set and function set
    nodes = Array{Float64}(undef, (0,0))
    f = Array{Float64}(undef, 0)

    # Initialize differentiation matrices
    Dx = Array{Float64}(undef, (0,0))
    Dy = Array{Float64}(undef, (0,0))
    
    # Initialize reinitialization counter
    reset = initCnt
    
    # Begin evolution in time
    for i∈time
        # Check for reinitialization
        if reset >= initCnt
            # Add node bands to node set and orient the surface
            (nodes,f) = generateNodeBand(zeroSet)
            if nodes[1,circ1[1]] - nodes[1,circ1[N+1]] < 0
                f[circ1] *= -1
            end
            # if nodes[circ2[floor(Int,(N+1)/2)]] - nodes[1,circ2[N+floor(Int,(N+1)/2)]] > 0
            #     f[circ2] *= -1
            # end
            if nodes[1,circ2[1]] - nodes[1,circ2[N+1]] < 0
                f[circ2] *= -1
            end
            
            # Discretize first order derivatives
            Dx = discretize∂xi(nodes, n, m, o, 1)
            Dy = discretize∂xi(nodes, n, m, o, 2)


            # Reset reinitialization counter
            reset = 0
        end

        # Compute size of node distribution
        NN = size(nodes,2)

        # Compute needed derivatives
        ϕx = Dx*f
        ϕy = Dy*f

        # Compute ∇ϕ
        ∇ϕ = [ϕx';ϕy']

        # Compute the normal direction: ∇ϕ/||∇ϕ||
        norm∇ϕ = zeros(2,NN)
        for j∈1:NN
            norm∇ϕ[:,j] = ∇ϕ[:,j]/norm(∇ϕ[:,j])
        end

        # Compute curvature or the divergence of the normal direction
        κ = abs.(Dx*norm∇ϕ[1,:] + Dy*norm∇ϕ[2,:])

        # Compute (normal) ⋅ ∇ϕ
        nDot∇ϕ = zeros(NN)
        for j∈1:NN
            nDot∇ϕ[j] = norm∇ϕ[:,j]⋅∇ϕ[:,j]
        end
        
        # Compute next time step
        f -= κ .* nDot∇ϕ * Δt

        # Find new zero level-set
        zeroSet = nodes[:,1:2N]
        zeroSet = newtonSolve(zeroSet, nodes, f, n=5, m=m, o=-1,
                              ε=10^(-13), maxIts=100)
        
        # Plot zero level-set
        plotA = scatter(nodes[1,:],nodes[2,:],
                        marker_z = f,
                        ratio = 1,
                        markersize = 2,
                        xlims = (-0.6,0.6),
                        ylims = (-0.6,0.6))
        # scatter!(zeroSet[1,:],zeroSet[2,:],
        #          ratio = 1,
        #          markersize = 1,
        #          xlims = (-0.6,0.6),
        #          ylims = (-0.6,0.6))
        display(plotA)
        sleep(0.01)

        # Increment reinitialization counter
        reset += 1
    end
end

grow(25, n=9, m=5, o=2, Δt=10^(-5), tf=0.01, initCnt = 200) 
