using Plots
using Analysis
using ColorSchemes
using LaTeXStrings
using SparseArrays
include("LevelSet.jl")

# Surface parameterization
x(r,t) = r*cos.(t)
y(r,t) = r*sin.(t)

# Radius over time
rt(r0,t) = sqrt(r0^2 + 2t)

function initLevelSet(surfNodes, backNodes, norms; m=5)
    N  = size(backNodes, 2)             # Number of background nodes
    sN = size(surfNodes, 2)             # Number of surface nodes
    D  = size(backNodes, 1)             # Number of dimensions/variables

    A     =  zeros(sN + D*N, N)      # Full constraint matrix
    CZero =  zeros(sN, N)            # Zero level set constraints
    CDxi  = [zeros(N,N) Pkgfor i ∈ 1:D] # Normal derivative constrains

    # Populate each constraint matrix
    for i ∈ 1:sN
        for j ∈ 1:N
            CZero[i, j] = ϕ(surfNodes[:, i], backNodes[:, j], m)
        end
    end
    for i ∈ 1:N
        for j ∈ 1:N
            for k ∈ 1:D
                CDxi[k][i, j] = ϕ_xi(backNodes[:, i], backNodes[:, j], m, k)
            end
        end
    end

    A[1:sN, :] = CZero
    for i = 0:D - 1
        A[(i * N + 1 + sN):(sN + (i + 1) * N), :] = CDxi[i + 1]
    end

    # Populate constraint vector
    b = zeros(D*N + sN)
    for i = 0:D - 1
        b[(i * N + 1 + sN):((i + 1) * N + sN)] = norms[i + 1,:]
    end

    λ = A \ b       # Compute pieceinterpolation weights

    f = zeros(N)    # Distance function
    for i = 1:N
        for j = 1:N
            f[i] += λ[j] * ϕ(backNodes[:, i], backNodes[:,j], m)
        end
    end

    return f
end

function initTest(sN,N; n, m, o, r, c, width)
    # Coordinate discretization
    t = range(0,2*π-2*π/sN, length = sN)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c; y(r,t)'.+c y(r,t)'.-c]

    # Generate node band
    nodes = hexBand(zeroSet, N=N, n=width)
    N = size(nodes,2)

    # Normal directions to the zero-set
    norms = zeros(2, N)
    for i ∈ 1:N
        tmp = [0,0]
        if nodes[2,i] >= -nodes[1,i]
            tmp = nodes[:, i] - [c,c]
        else
            tmp = nodes[:, i] + [c,c]
        end
        norms[:, i] = tmp / norm(tmp)
    end

    # Initialize distance function
    f = initLevelSet(zeroSet, nodes, norms, m=m)

    plotA = scatter(nodes[1,:], nodes[2,:],
                    marker_z    = f,
                    c           = :viridis,
                    legend      = false,
                    colorbar    = true,
                    markeralpha = 0.75,
                    markersize  = 3,
                    ratio       = 1)
    display(plotA)
end

function timeSolve(sN,N; n, m, o, r, c, width, Δt, tf)
    # Coordinate discretization
    t = range(0,2*π-2*π/sN, length = sN)

    # Compute surface nodes
    zeroSet = [x(r,t)'.+c x(r,t)'.-c; y(r,t)'.+c y(r,t)'.-c]

    # Generate node band
    nodes = hexBand(zeroSet, N=N, n=width)
    N = size(nodes,2)

    # Normal directions to the zero-set
    norms = zeros(2, N)
    for i ∈ 1:N
        tmp = [0,0]
        if nodes[2,i] >= -nodes[1,i]
            tmp = nodes[:, i] - [c,c]
        else
            tmp = nodes[:, i] + [c,c]
        end
        norms[:, i] = tmp / norm(tmp)
    end

    # Initialize distance function
    f = initLevelSet(zeroSet, nodes, norms, m=m)

    plotA = scatter(nodes[1,:], nodes[2,:],
                    marker_z    = f,
                    c           = :viridis,
                    legend      = false,
                    colorbar    = true,
                    markeralpha = 0.75,
                    markersize  = 3,
                    ratio       = 1)
    display(plotA)

    # Descretize required operators
    Dx = discretize∂xi(nodes, n, m, o, 1)
    Dy = discretize∂xi(nodes, n, m, o, 2)

    # Define RHS operator
    function RHS!(du, u)
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

    # Time variable
    t = [0:Δt:tf...]

    # Solve in time
    df = zeros(N)
    for i ∈ 2:length(t)
        # Step in time
        RHS!(df, f)
        f += df * Δt

        # Check for plotting and reinitialization
        if i % 50 == 0
            # Compute zeroSet for solution visualization and analysis
            zeroSet[:, :] = coulNewtonAdapt(zeroSet, nodes, f,
                                      n=n, m=m, o=o,
                                      maxIts=500, μ=1.25,
                                      Δt=10^(-5), ε=10^(-13))

            # Reinitialize nodes
            fx = Dx*f
            fy = Dy*f
            for j = 1:N
                norms[:,j] = [fx[j],fy[j]]
                norms[:,j] /= norm(norms[:,j])
            end
            f = initLevelSet(zeroSet, nodes, norms, m = m)

            # Compute true radius of circles at current time
            rts = fill(rt(r, t[i]), 2sN)

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
        end
    end
end

# initTest(50,2000, n=15, m=5, o=2, r=1, c=0.8, width=60)
timeSolve(100,3000, n=15, m=5, o=2, r=1, c=0.8, width=80, Δt = 10^(-4), tf = 0.05)
