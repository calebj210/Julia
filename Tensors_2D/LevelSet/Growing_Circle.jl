### Test for solving the growing circle in time

# Packages to include and use
include("LevelSet.jl")
using Plots

# Smoothed sign function
sgn(ϕ,ε) = ϕ/sqrt(ϕ^2+ε^2)

function grow1(N,M; n, m, o, r, Δt, tf)
    # Form our grid
    x = -1.25r:2.5r/(N-1):1.25r
    y = -1.25r:2.5r/(M-1):1.25r
    nodes = [repeat(x,inner = (M,1))';
             repeat(y,outer = (N,1))']

    # Compute grid spacing
    ε = 10(x[2]-x[1])

    # Level set surface
    f = zeros(N*M)
    for i∈1:N*M
        f[i] = norm(nodes[:,i]) - r
    end

    # Plot our initial surface
    plotA = plot(x, y, reshape(f,N,M),
                 st = :contour,
                 levels = [0,-0.25,0.25,0.5],
                 ratio = 1)

    display(plotA)

    # Descretize required operators
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)
    Dxx = discretize∂xii(nodes, n, m, o, 1)
    Dyy = discretize∂xii(nodes, n, m, o, 2)
    Dxy = discretize∂x∂y(nodes, n, m, o)

    # for i∈1:M
    #     Dx[i,:] .= 0
    #     Dx[i,i] = 1
    #     Dy[i,:] .= 0
    #     Dy[i,i] = 1
    #     Dxx[i,:] .= 0
    #     Dxx[i,i] = 1
    #     Dyy[i,:] .= 0
    #     Dyy[i,i] = 1
    #     Dxy[i,:] .= 0
    #     Dxy[i,i] = 1
    # end
    # for i∈M:M:N*M
    #     Dx[i,:] .= 0
    #     Dx[i,i] = 1
    #     Dy[i,:] .= 0
    #     Dy[i,i] = 1
    #     Dxx[i,:] .= 0
    #     Dxx[i,i] = 1
    #     Dyy[i,:] .= 0
    #     Dyy[i,i] = 1
    #     Dxy[i,:] .= 0
    #     Dxy[i,i] = 1
    # end
    # for i∈(N-1)*M+1:N*M
    #     Dx[i,:] .= 0
    #     Dx[i,i] = 1
    #     Dy[i,:] .= 0
    #     Dy[i,i] = 1
    #     Dxx[i,:] .= 0
    #     Dxx[i,i] = 1
    #     Dyy[i,:] .= 0
    #     Dyy[i,i] = 1
    #     Dxy[i,:] .= 0
    #     Dxy[i,i] = 1
    # end
    # for i∈1:M:N*M-1
    #     Dx[i,:] .= 0
    #     Dx[i,i] = 1
    #     Dy[i,:] .= 0
    #     Dy[i,i] = 1
    #     Dxx[i,:] .= 0
    #     Dxx[i,i] = 1
    #     Dyy[i,:] .= 0
    #     Dyy[i,i] = 1
    #     Dxy[i,:] .= 0
    #     Dxy[i,i] = 1
    # end
    
    # Initialize time variable
    time = Δt:Δt:tf

    # Initialize reinitialization counter
    reset = 0
    
    # Begin evolution in time
    for i∈time
        # # Check for reinitialization
        # if reset == 1
        #     for j∈1:100
        #         # Compute ∂x and ∂y
        #         dx = Dx*f
        #         dy = Dy*f
                
        #         # Compute gradient magnitudes
        #         mag = zeros(N*M)
        #         for k∈1:N*M
        #             mag[k] = norm([dx[k],dy[k]])
        #         end
            
        #         # Solve gradient scaling equation
        #         f -= sgn.(f,ε).*(mag .- 1)*10^(-4)
        #     end
        # else
        #     reset += 1
        # end
        
        # Compute need derivatives
        ϕx = Dx*f
        ϕy = Dy*f
        ϕxx = Dxx*f
        ϕyy = Dyy*f
        ϕxy = Dxy*f

        # Compute curvature
        κ = -abs.(ϕxx.*ϕy.^2 - 2ϕx.*ϕy.*ϕxy + ϕyy.*ϕx.^2)./(ϕx.^2 +
                                                           ϕy.^2)
        
        # Compute next time step
        f += κ*Δt

        plotA = plot(x, y, reshape(f,N,M),
                     st = :contour,
                     levels = [0],
                     ratio = 1)
        display(plotA)
    end
end

function grow2(N,M; n, m, o, r, Δt, tf)
    # Form our grid
    x = -2r:4r/(N-1):2r
    y = -2r:4r/(M-1):2r
    nodes = [repeat(x,inner = (M,1))';
             repeat(y,outer = (N,1))']

    # Level set surface
    f = zeros(N*M)
    for i∈1:N*M
        f[i] = norm(nodes[:,i]) - r
    end

    # Plot our initial surface
    plotA = plot(x, y, reshape(f,N,M),
                 st = :contour,
                 levels = [0,-0.25,0.25,0.5],
                 ratio = 1)
    display(plotA)

    # Descretize required operators
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)
    
    # Initialize time variable
    time = Δt:Δt:tf

    # Initialize reinitialization counter
    reset = 0
    
    # Begin evolution in time
    for i∈time
        # Compute needed derivatives
        ϕx = Dx*f
        ϕy = Dy*f

        # Compute ∇ϕ
        ∇ϕ = [ϕx';ϕy']

        # Compute the normal direction: ∇ϕ/||∇ϕ||
        norm∇ϕ = zeros(2,N*M)
        for j∈1:N*M
            norm∇ϕ[:,j] = ∇ϕ[:,j]/norm(∇ϕ[:,j])
        end

        # Compute curvature or the divergence of the normal direction
        κ = abs.(Dx*norm∇ϕ[1,:] + Dy*norm∇ϕ[2,:])

        # Compute normal ⋅ ∇ϕ
        nDot∇ϕ = zeros(N*M)
        for j∈1:N*M
            nDot∇ϕ[j] = norm∇ϕ[:,j]⋅∇ϕ[:,j]
        end
        
        # Compute next time step
        f -= κ .* nDot∇ϕ * Δt

        plotA = plot(x, y, reshape(f,N,M),
                     st = :contour,
                     levels = [0],
                     ratio = 1)
        display(plotA)
    end
end

## Circle parameterization
x(t,r) = r*cos.(t)
y(t,r) = r*sin.(t)

function grow3(N; n, m, o, r, Δt, tf, initCnt)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Creating zero level-set
    zeroSet = [x(t,r)';y(t,r)']
    
    # Plot full node set
    plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                     ratio = 1,
                     markersize = 2,
                    xlims = (-2r,2r),
                    ylims = (-2r,2r))
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
            if nodes[1,1] - nodes[1,N+1] < 0
                f *= -1
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
        zeroSet = nodes[:,1:N]
        zeroSet = newtonSolve(zeroSet, nodes, f, n=5, m=m, o=-1,
                              ε=10^(-13), maxIts=5)
        
        # Plot zero level-set
        # plotA = scatter(nodes[1,:],nodes[2,:],
        #                 ratio = 1,
        #                 markersize = 1,
        #                 xlims = (-2r,2r),
        #                 ylims = (-2r,2r))
        # scatter!(zeroSet[1,:],zeroSet[2,:],
        #          ratio = 1,
        #          markersize = 1,
        #          xlims = (-2r,2r),
        #          ylims = (-2r,2r))
        # display(plotA)

        # Increment reinitialization counter
        reset += 1
    end
end

# grow1(10,10, n=100, m=5, o=2, r=0.2, Δt=10^(-5), tf=0.02)

# grow2(10,10, n=23, m=5, o=2, r=0.2, Δt=10^(-5), tf=0.02)

grow3(35, n=9, m=5, o=2, r=0.2, Δt=10^(-5), tf=0.06, initCnt = 25) 
