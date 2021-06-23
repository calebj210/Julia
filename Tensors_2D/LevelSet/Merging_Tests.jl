### Merging Circles Test

# Needed packages
using Plots
using Analysis
using ColorSchemes
using LaTeXStrings
using DifferentialEquations
include("LevelSet.jl")

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Radius over time
rt(r0,t) = sqrt(r0^2 + 2t)

# Test with hexagonal grid and no reinitialization
function simulate(sN,N; n, m, o, r, c, maxIts=250, Δt, tf, width)
    # Coordinate discretization
    t = range(0,2*π-2*π/sN, length = sN)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c;y(r,t)'.+c y(r,t)'.-c]

    # Generate node band
    nodes = hexBand(zeroSet, N=N, n=width)
    N = size(nodes,2)

    # Plot hex band
    plotA = scatter(nodes[1,:],nodes[2,:],
                    legend = false,
                    colorbar = true,
                    markeralpha = 0.75,
                    markersize = 2,
                    ratio = 1,
                    dpi = 300)
    display(plotA)

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

    ### Uncomment to show eigenvalues and condition numbers of Dx and Dy
    # show(eigvals(Array(Dx)))
    # show(eigvals(Array(Dy)))
    # print("cond(Dx) = ", log10(cond(Array(Dx))), ".\n")
    # print("cond(Dy) = ", log10(cond(Array(Dy))), ".\n")

    # Initialize time variable
    time = Δt:Δt:tf

    # Initialize reinitialization counter and error storage
    k = 1
    maxErr = 0
    errs = zeros(sN)

    # Plot initial surface
    # plotA = scatter(zeroSet[1,:],zeroSet[2,:],
    #                 ratio = 1,
    #                 legend = false,
    #                 xlims = (-(c+2r),c+2r),
    #                 ylims = (-(c+2r),c+2r),
    #                 markersize = 2,
    #                 dpi=300)
    # display(plotA)

    # Begin time evolution
    for i∈time
        # Compute needed derivatives
        ϕx = Dx*f
        ϕy = Dy*f

        # Compute laplacian hyperviscosity
        # Δϕ = 0*f
        Δϕ = 10^(-2)*(Dx*Dx + Dy*Dy)^2*f

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
        f -= (κ .* nDot∇ϕ) * Δt

        # Check for solution plotting and error calulations
        if i == time[end]
            # Compute zeroSet for solution visualization and analysis
            zeroSet = coulNewtonAdapt(zeroSet, nodes, f, n=n, m=m, o=o, maxIts=500, μ=1.25, Δt=10^(-5), ε=10^(-13))

            # # Plot solution
            plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                            ratio = 1,
                            legend = false,
                            xlims = (-(c+2r),c+2r),
                            ylims = (-(c+2r),c+2r),
                            markersize = 2,
                            dpi=300)

            # Plot hex band
            # plotA = scatter(nodes[1,:],nodes[2,:],
            #         legend = false,
            #         marker_z = f,
            #         colorbar = true,
            #         markeralpha = 0.75,
            #         markersize = 4,
            #         ratio = 1,
            #         dpi = 300)
            #         display(plotA)

            # Compute true radius of circles at current time
            rts = fill(rt(r, i), 2sN)

            # Compute radius of each node
            rads = zeros(2sN)
            for j∈1:2sN
                if zeroSet[2,j] >= -zeroSet[1,j]
                    rads[j] = norm(zeroSet[:,j] - [c,c])
                else
                    rads[j] = norm(zeroSet[:,j] + [c,c])
                end
            end

            # Compute errors
            maxErr, errs = scalError(rts,rads)

            # Plot solution with errors overlayed
            display(scatter(zeroSet[1,:],zeroSet[2,:],
                            marker_z = log10.(errs),
                            c = :viridis,
                            legend = false,
                            colorbar = true,
                            markeralpha = 0.75,
                            markersize = 3,
                            ratio = 1))
        end

        k += 1
    end

    # Print errors to terminal
    display(norm(errs))

    # return (maxErr, errs)
    return [N, maxErr]
end

# Node density tests
function densityTest(N, n)
    # Number of tests
    tN = size(N,1)

    # Preallocate space for errors
    errs = Array{Any}(undef,(3,tN))

    # Run tests and record errors
    for i∈eachindex(N)
        # Compute width
        wdth = ceil(Int, n*N[i]/N[1])

        # Compute error with current density
        errs[1:2,i] = simulate(100, N[i], n=15, m=3, o=2, r=1, c=0.8, Δt=10^(-5), tf=0.05, width=wdth)

        # Did the circles merge
        print("Did the circles merge?\n")
        errs[3,i] = readline()
    end

    return errs
end

# Band width test
function widthTest(N, n)
    # Number of tests
    tN = size(n,1)

    # Preallocate space for erros
    errs = Array{Any}(undef, (3,tN))

    # Run tests and record errors
    for i∈eachindex(n)
        # Compute errors
        errs[1:2,i] = [n[i], simulate(100, N, n=15, m=3, o=2, r=1, c=0.8, Δt=10^(-5), tf=0.05, width=n[i])[2]]

        # Did the circles merge
        print("Did the circles merge?\n")
        errs[3,i] = readline()
    end

    return errs
end

# Time step test
function timeTest(N, n, Δt)
    # Number of tests
    tN = size(Δt,1)

    # Preallocate space for erros
    errs = Array{Any}(undef, (3,tN))

    # Run tests and record errors
    for i∈eachindex(Δt)
        # Compute errors
        errs[1:2,i] = [Δt[i], simulate(100, N, n=15, m=3, o=2, r=1, c=0.8, Δt=Δt[i], tf=0.05, width=n)[2]]

        # Did the circles merge
        print("Did the circles merge?\n")
        errs[3,i] = readline()
    end

    return errs
end

# Run single test
simulate(100,2000, n=21, m=5, o=2, r=1, c=0.8, Δt=10^(-4), tf=0.05, width=10)

# Run density test
# test1 = densityTest([250,500,1000,2000,4000],5)

# Run band width test
# test2 = widthTest(1000, [3,6,12,24,48,96])

# Run time step test
# test3 = timeTest(1000, 15, 0.5 .^ [0,1,2,3,4,5] * 0.25*10^(-3))
