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
function simulate(sN,N; n, m, o, r, c, maxIts=250, tf, width)
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

    # Descretize required operators
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)

    # Define RHS operator
    function RHS!(du, u, p, t)
        # Compute needed derivatives
        ϕx = Dx*u
        ϕy = Dy*u

        # Compute ∇ϕ / gradient of u
        ∇ϕ = [ϕx';ϕy']

        # Compute the normal direction: ∇ϕ/||∇ϕ||
        norm∇ϕ = zeros(2,N)
        for j∈1:N
            norm∇ϕ[:,j] = ∇ϕ[:,j]/norm(∇ϕ[:,j])
        end

        # Compute curvature / the divergence of the normal direction
        κ = abs.(Dx*norm∇ϕ[1,:] + Dy*norm∇ϕ[2,:])

        # Compute normal ⋅ ∇ϕ
        nDot∇ϕ = zeros(N)
        for j∈1:N
            nDot∇ϕ[j] = norm∇ϕ[:,j]⋅∇ϕ[:,j]
        end

        # Compute RHS
        du[1:end] = -κ .* nDot∇ϕ
    end

    # Initialize node band
    f = zeros(N)
    for i ∈ 1:N
        if nodes[2,i] >= -nodes[1,i]
            f[i] = norm(nodes[:,i] - [c,c]) - r
        else
            f[i] = norm(nodes[:,i] + [c,c]) - r
        end
    end

    # Initialize reinitialization counter and error storage
    maxErr = 0
    errs = zeros(sN)

    # Solve ODE in time
    prob = ODEProblem(RHS!, f, (0.0, tf))
    sol = solve(prob)


    # Plot solution
    for t ∈ 0.0:tf/4:tf
        # Compute zeroSet for solution visualization and analysis
        zeroSet = coulNewtonAdapt(zeroSet, nodes, sol(t),
                                  n=n, m=m, o=o,
                                  maxIts=500, μ=1.25,
                                  Δt=10^(-5), ε=10^(-13))

        # Compute true radius of circles at current time
        rts = fill(rt(r, t), 2sN)

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
        plotA = scatter(zeroSet[1,:],zeroSet[2,:],
                        marker_z = log10.(errs),
                        c = :viridis,
                        legend = false,
                        colorbar = true,
                        markeralpha = 0.75,
                        markersize = 3,
                        ratio = 1)

        display(plotA)

        sleep(0.5)
    end
end

# Run single test
simulate(50,2000, n=15, m=5, o=2, r=1, c=0.8, tf=0.05, width=60)
