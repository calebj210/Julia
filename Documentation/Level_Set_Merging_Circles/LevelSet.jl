### Level Set Functions
# Including used packages
using SparseArrays
include("LevelSetup.jl")
include("RBF_Functions.jl")

# Function for generating node band around surface
function generateNodeBand(Nodes)
    # Compute the number of nodes
    N = size(Nodes,2)

    # Find nearest neighbors
    idx = knnFull(Nodes, 5)

    # Find the approximate normals of our surface
    norms = approxNormals(Nodes, idx)

    # Orient surface normals
    norms = orVecs(Nodes, norms, idx)

    # Scale normals for adding nodes
    for i ∈ 1:N
        norms[:,i] *= 0.25*norm(Nodes[:,i] - Nodes[:,idx[i][2]])
    end

    # Compute new node set with bands
    nodes = [Nodes (Nodes - norms) (Nodes + norms) (Nodes - 2norms) (Nodes + 2norms)]

    # Assign distance values to the node set
    # f = [fill(0,N); fill(-1,N); fill(1,N)]
    f = zeros(5N)
    for i∈1:N
        tmp = norm(norms[:,i])
        f[N+i] = -tmp
        f[2N+i] = tmp
        f[3N+i] = -2tmp
        f[4N+i] = 2tmp
    end

    return (nodes,f)
end

# Hexagonal grid generator
function hexGen(N; minx, maxx, miny, maxy)
    # Compute bound ranges
    Δx = maxx - minx
    Δy = maxy - miny

    # Compute number of nodes in each dimension
    n = round(Int,sqrt(N*Δx/Δy))
    m = round(Int, sqrt(N*Δy/Δx))
    
    # Compute base grid
    x = 0:n-1
    y = 0:m-1
    nodes = [repeat(x, inner = (m,1))' / ((n-1)+0.5) * Δx;
             repeat(y, outer = (n,1))' / (m-1)*Δy]

    # Compute x offsets
    xOff = repeat(repeat([0,0.5 / ((n-1)+0.5) * Δx],
                         outer = (ceil(Int,m/2),1))[1:m],
                  outer = (n,1)) .+ minx

    # Offset the nodes
    nodes[1,:] += xOff
    nodes[2,:] .+= miny
    
    return nodes
end

# Hexagonal band generartor
function hexBand(zeroSet; N=1000, n=15)
    # Compute node bounds
    mins = minimum(zeroSet, dims=2)
    maxs = maximum(zeroSet, dims=2)
    
    # Generate background hex grid
    nodes = hexGen(N, minx=2mins[1], maxx=2maxs[1], miny=2mins[2], maxy=2maxs[2])

    # Generate KDTree of node set
    kdtree = KDTree(nodes)

    # Compute nearest neighbors
    idx = knn(kdtree, zeroSet, n, false)[1]

    # Take union of nearest neighbors
    idx = union(idx...)

    # Remove nodes not near the surface
    nodes = nodes[:,idx]

    return nodes
end
function hexBand(zeroSet, bNodes; n=15)
    # Generate KDTree of node set
    kdtree = KDTree(bNodes)

    # Compute nearest neighbors
    idx = knn(kdtree, zeroSet, n, false)[1]

    # Take union of nearest neighbors
    idx = union(idx...)

    # Remove nodes not near the surface
    nodes = bNodes[:,idx]
end

# Distance function initializer
function initialize(surfNodes, backNodes, norms; n=5, m=5, o=0, maxIts=500, ε=10^(-10), Δt=1, prompt=false)
    # Compute number of nodes to be initialized
    N = size(backNodes, 2)

    # Compute KDTree for surface node set
    surfKDTree = KDTree(surfNodes)
    
    # Compute nearest neighbors among the surface nodes
    surfIdx = knn(surfKDTree, surfNodes, n, true)[1]

    # Compute the nearest surface node frome each node in the background node set     
    backIdx = nn(surfKDTree, backNodes)[1]

    # Construct polynomial matrix
    polMat = polyMat(o,1)
    
    # Initialize distance function
    F = zeros(N)

    # Compute distance from each node in the background to the surface
    for i∈1:N
        # Center surface nodes and node of interest at origin about nearest surface node
        nodes = cent([surfNodes[:,surfIdx[backIdx[i]]] backNodes[:,i]])
        
        # Rotate centered nodes about the origin to face "up" in the normal direction
        rotUp!(nodes, norms[:,backIdx[i]])

        # Interpolate over the nodes
        λ = findλ(nodes[1,1:n]',nodes[2,1:n], polMat, m=m)

        # Select center node as initial guess for closest approach distance
        x = 0.0

        # Initialize surface function value and derivative
        f = nodes[2,1]
        f_x = 0.0

        # Initialize criteria
        xStep = Inf

        # Begin Newton's minimization algorithm
        for j∈1:maxIts
            # Compute surface function and its derivatives at current x
            f = S(nodes[1,1:n]', [x], λ, polMat, m=m)
            f_x = S_xi(nodes[1,1:n]', [x], λ, 1, polMat, m=m)
            f_xx = S_xii(nodes[1,1:n]', [x], λ, 1, polMat, m=m)
            
            # Exit loop if convergence criteria is met
            if abs(xStep) < ε
                if prompt
                    display(string("Converged in ", j, " iterations"))
                end
                break
            end

            # Compute distance function derivatives
            d_x = -2(nodes[1,n+1]-x + f_x*(nodes[2,n+1]-f))
            d_xx = 2(1-f_xx*(nodes[2,n+1]-f) + f_x^2)

            # Compute x-step
            xStep = d_x/d_xx
            
            # Iterate
            x -= Δt*xStep
        end

        # Compute magnitude of the distance function
        F[i] = norm(nodes[:,n+1] - [x,f])
        
        # Compute the scaled normal direction
        nrml = [-f_x, 1]
        
        # Sign the distance function
        F[i] *= sign(nrml⋅(nodes[:,n+1]-nodes[:,1]))
    end

    return F
end

# Node reinitialization algorithm
function reinit(nodes; oldNodes, F, zeroSet, smooth = false, n=1, m=5, o=0)
    # Number of nodes
    N = size(nodes, 2)

    # Generate KD trees of node sets
    zeroKDTree = KDTree(zeroSet)
    oldKDTree = KDTree(oldNodes)

    # Compute nearest neighbors to the newest node set
    dists = nn(zeroKDTree, nodes)[2]
    idx = knn(oldKDTree, nodes, n, true)[1]

    # Compute distance function
    f = zeros(N)
    for i∈1:N
        f[i] = sign(F[idx[i][1]])*dists[i][1]
    end

    # Check for node replacement
    if smooth
        # Compute nearest neighbors to the zero-set
        newKDTree = KDTree(nodes)
        zeroIdx = knn(newKDTree, zeroSet, 3, true)[1]

        # Construct the polynomial matrix
        polMat = polyMat(o)

        # Compute new function values at each node around the zero-set
        for i∈eachindex(zeroSet[1,:]), j∈1:3
            # Compute current index
            k = zeroIdx[i][j]
            
            # Compute interpolant weights
            λ = findλ(oldNodes[:,idx[k]], F[idx[k]], polMat, m=m)

            # Compute new function value
            f[k] = S(oldNodes[:,idx[k]], nodes[:,k], λ, polMat, m=m)
        end
    end

    return f
end

# Coul-Newton method for adaptive nodes
function coulNewtonAdapt(zrs, Nodes, Finit; n, m, o, maxIts=200, μ=2, Δt=0.1, ε=10^(-10))
    # Number of nodes
    N = size(Nodes,2)
    # Number of zeros
    Z = size(zrs,2)
    
    # Compute polynomial terms and number of terms
    polMat = polyMat(o)
    P = size(polMat,2)

    # Initialize f
    f = [Finit;zeros(N+Z)]

    # Combine zero-level-set and nodes
    nodes = [Nodes zrs]
    
    # Compute KDTree for background nodes
    kdtree = KDTree(nodes[:,1:N])
    
    # Begin Coulomb-Newton iteration on zero-level-set
    for i∈1:maxIts
        # display(scatter(nodes[1,N+1:end],nodes[2,N+1:end],
        #                 ratio = 1))
        # sleep(.01)

        # Compute nearest neigbors
        idx, tmp = knn(kdtree, nodes[:,N+1:end], n, true)

        # Compute nearest neighbors on the zero-level-set
        zeroKDTree = KDTree(nodes[:,N+1:end])
        zerIdx, dists = knn(zeroKDTree, nodes[:,N+1:end], 5, true)
        
        # Reset Coulombic force and f*∇f/||∇f||
        Fc = zeros(2,Z)
        ∇f = zeros(2,Z)
        
        # Compute gradient and function values at each node as well as the
        # Coulombic forces
        for j∈1:Z
            # Compute interpolant weights
            λ = findλ(nodes[:,idx[j]], f[idx[j]], polMat, m=m)
            
            # Compute function values and first order derivatives
            f[N+j] = S(nodes[:,idx[j]], nodes[:,N+j], λ, polMat, m=m)
            fx = S_xi(nodes[:,idx[j]], nodes[:,N+j], λ, 1, polMat, m=m)
            fy = S_xi(nodes[:,idx[j]], nodes[:,N+j], λ, 2, polMat, m=m)

            # Compute f*∇f/||∇f||
            ∇f[:,j] = f[N+j]*[fx,fy]/(norm([fx,fy])^2)

            # Compute Coulombic force
            for k∈2:5
                Fc[:,j] += (nodes[:,N+j] - nodes[:, N+zerIdx[j][k]]) / dists[j][k]^3
            end
        end

        # Scale Coulombic force gradient
        Fc *= μ*Δt/2

        # Iterate nodes
        nodes[:,N+1:end] += Fc - ∇f
    end

    zrs = newtonSolve(nodes[:,N+1:end], nodes[:,1:N], f[1:N], n=n, m=m, o=o,
    maxIts = 200, ε = ε) 
    
    return zrs
end


# Coulomb-Newton Method for generating node band
function coulNewtonBand(Nodes, Finit; n, m, o, maxIts=200, μ=2, η=0.01, Δt=0.1, ε=10^(-3))
    # Number of nodes
    N = size(Nodes,2)
    
    # Compute polynomial terms and number of terms
    polMat = polyMat(o)
    P = size(polMat,2)

    # Initialize f
    f = copy(Finit)

    # Combine zero-level-set and nodes
    nodes = copy(Nodes)
    
    # Compute KDTree for background nodes
    kdtree = KDTree(nodes)

    # Compute ambient node bounds and ranges
    xmin = minimum(nodes[1,:])
    xmax = maximum(nodes[1,:])
    xrng = xmax-xmin
    
    ymin = minimum(nodes[2,:])
    ymax = maximum(nodes[2,:])
    yrng = ymax-ymin
    
    # Construct initial node band
    xmin += xrng/8
    xmax -= xrng/8
    ymin += yrng/8
    ymax -= yrng/8
    
    band = hexGen(N, minx=xmin, maxx=xmax, miny=ymin, maxy=ymax)

    # Put random noise into the node band
    tmp = 0.25norm(band[:,1]-band[:,2])
    for i∈eachindex(band[1,:])
        band[:,i] += 2tmp*(rand(2) .- 0.5)
    end

    # Compute number of nodes in band
    bN = size(band,2)

    # plotA = scatter(nodes[1,:],nodes[2,:],
    #                 marker_z = Finit,
    #                 c = :viridis,
    #                 legend = false,
    #                 colorbar = true,
    #                 colorbar_title = "Distance Function",
    #                 markeralpha = 0.9,
    #                 markerstrokewidth = 0,
    #                 ratio = 1,
    #                 markersize = 3,
    #                 markershape = :+,
    #                 xlims = (-2,2),
    #                 ylims = (-2,2),
    #                 xticks = [-2:0.5:2...],
    #                 yticks = [-2:0.5:2...])
    # scatter!(band[1,:],band[2,:],
    #         c = :grey,
    #          markeralpha = 0.5,
    #          markersize = 1.5)
    # display(plotA))

    # Preallocate space for node band function values
    bf = Array{Float64}(undef, bN)
    
    # Begin Coulomb-Newton iteration on background nodes
    for i∈1:maxIts
        # Compute KDTree for nodes
        kdtree = KDTree(nodes)

        # Compute KDTree for node band
        bandKDTree = KDTree(band)
        
        # Compute nearest neighbors of node band to the ambient node set
        idx, tmp = knn(kdtree, band, n, true)

        # Compute nearest neighbors and their distances in the node band
        bandIdx, dists = knn(bandKDTree, band, n, true)
        
        # Reset Coulombic force and f*∇f/||∇f||
        Fc = zeros(2,bN)
        ∇f = zeros(2,bN)
        
        # Compute gradient and function values at each node as well as the
        # Coulombic forces
        for j∈1:bN
            # Compute interpolant weights
            λ = findλ(nodes[:,idx[j]], f[idx[j]], polMat, m=m)
            
            # Compute function values and first order derivatives
            bf[j] = S(nodes[:,idx[j]], band[:,j], λ, polMat, m=m)
            fx = S_xi(nodes[:,idx[j]], band[:,j], λ, 1, polMat, m=m)
            fy = S_xi(nodes[:,idx[j]], band[:,j], λ, 2, polMat, m=m)

            # Compute f*∇f/||∇f||
            ∇f[:,j] = bf[j]*[fx,fy]/(norm([fx,fy])^2)

            # Compute Coulombic force
            for k∈2:n
                Fc[:,j] += (band[:,j] - band[:, bandIdx[j][k]]) / dists[j][k]^2
            end
        end

        # Scale Coulombic force gradient
        Fc *= μ*Δt
        ∇f *= η

        # Iterate nodes
        band += Fc - ∇f
        
        # plotA = scatter(band[1,:],band[2,:],
        #                 c = :grey,
        #                 markeralpha = 0.75,
        #                 markerstrokewidth = 0,
        #                 ratio = 1,
        #                 markersize = 1.5,
        #                 xlims = (-1.25,1.25),
        #                 ylims = (-1.25,1.25),
        #                 dpi = 300)
        # scatter!(nodes[1,:],nodes[2,:],
        #          marker_z = Finit,
        #          c = :viridis,
        #          legend = false,
        #          colorbar = true,
        #          colorbar_title = "Distance Function",
        #          markeralpha = 0.5,
        #          markerstrokewidth = 0,
        #          ratio = 1,
        #          markersize = 3,
        #          markershape = :+)
        # display(plotA)

        # Check for convergence
        if maximum(abs.(Fc-∇f)/xrng) <= ε
            return (band, bf, bN)
        end
    end
    
    return (band, bf, bN)
end

# Routine for finding the normals
function findNormals(nodes, f, n, m, o)
    # Number of nodes
    N = size(nodes,2)

    # Compute polynomial degree matrix
    polMat = polyMat(o)

    # Compute nearest neighbors
    idx = knnFull(nodes, n)

    # Preallocate space for the normal directions
    norms = zeros(2,N)
    
    # Compute gradient of our level set at each node
    for i ∈ 1:N
        # Compute interpolant weights
        λ = findλ(nodes[:,idx[i]], f[idx[i]], polMat, m = m)

        # Compute gradient
        x1Grad = S_xi(nodes[:,idx[i]], nodes[:,i], λ, 1, polMat, m = m) 
        x2Grad = S_xi(nodes[:,idx[i]], nodes[:,i], λ, 2, polMat, m = m)
        grad = [x1Grad,x2Grad]

        # Normalize and store normals
        norms[:,i] = grad/norm(grad)
    end

    return norms
end

# Routine for computing the curvature
function computeκ(nodes, f, n, m, o)
    # Number of nodes
    N = size(nodes,2)

    # Compute polynomial degree matrix
    polMat = polyMat(o)

    # Compute nearest neighbors
    idx = knnFull(nodes, n)

    # Preallocate space for the normal directions
    κ = zeros(N)
    
    # Compute gradient of our level set at each node
    for i ∈ 1:N
        # Compute interpolant weights
        λ = findλ(nodes[:,idx[i]], f[idx[i]], polMat, m = m)

        # Compute all of the needed derivatives
        Dx1 = S_xi(nodes[:,idx[i]], nodes[:,i], λ, 1, polMat, m = m)
        Dx2 = S_xi(nodes[:,idx[i]], nodes[:,i], λ, 2, polMat, m = m)
        Dx11 = S_xii(nodes[:,idx[i]], nodes[:,i], λ, 1, polMat, m = m)
        Dx22 = S_xii(nodes[:,idx[i]], nodes[:,i], λ, 2, polMat, m = m)
        Dx12 = S_xij(nodes[:,idx[i]], nodes[:,i], λ, polMat, m = m)

        # Compute curvature
        κ[i] = (Dx11*Dx2^2 - 2Dx1*Dx2*Dx12 + Dx22*Dx1^2)/((Dx1^2 + Dx2^2)^(3/2))
    end

    return κ
end

# Routine for discretizing the Laplacian ∇²
function discretize∇²(nodes, n, m, o)
    # Number of nodes
    N = size(nodes,2)

    # Compute nearest neighbors
    idx = knnFull(nodes,n)

    # Compute polynomial degree matrix
    polMat = polyMat(o)
    
    # Number of polynomial terms
    P = size(polMat,2)

    # Preallocate space for the discretized operator matrix and the
    # linear operator vector
    D = zeros(N,N)
    b1 = zeros(n)
    b2 = zeros(P)
    
    # Begin populating discretized operator matrix
    for i∈1:N
        # Store node of interest
        X = nodes[:,i]

        # Compute collocation matrix
        A = colloc(nodes[:,idx[i]], polMat, m = m)
        
        ## Construct linear operator vector
        # Compute PHS components
        for j∈1:n
            b1[j] = ϕ_xii(X, nodes[:,idx[i][j]], m, 1) + ϕ_xii(X, nodes[:,idx[i][j]], m, 2)
        end

        # Compute polynomial components
        for j∈1:P
            # Store polynomial term degrees
            p = polMat[:,j]

            # Compute ∂²/∂x1² of polynomial
            if p[1]-2 >= 0
                tmpx1 = p[1]*(p[1]-1)X[1]^(p[1]-2)*X[2]^p[2]
            else
                tmpx1 = 0
            end

            # Compute ∂²/∂x2² of polynomial
            if p[2]-2 >= 0
                tmpx2 = p[2]*(p[2]-1)X[1]^p[1]*X[2]^(p[2]-2)
            else
                tmpx2 = 0
            end

            # Compute polynomial entry
            b2[j] = tmpx1 + tmpx2
        end

        # Combine PHS and polynomial of linear operator vector
        b = [b1;b2]
        
        # solve for the ith node weights
        λ = A\b

        # Store weights in the discretized operator matrix
        D[i,idx[i]] = λ[1:n]
    end

    return D
end

# ∂x discretizer 
function discretize∂xi(nodes, n, m, o, ii)
    # Number of nodes
    N = size(nodes,2)

    # Compute nearest neighbors
    idx = knnFull(nodes,n)

    # Compute polynomial degree matrix
    polMat = polyMat(o)
    
    # Number of polynomial terms
    P = size(polMat,2)

    # Preallocate space for the discretized operator matrix and the
    # linear operator vector
    Dx = spzeros(N,N)
    b1 = zeros(n)
    b2 = zeros(P)
    
    # Begin populating discretized operator matrix
    for i∈1:N
        # Store node of interest
        X = nodes[:,i]

        # Compute collocation matrix
        A = colloc(nodes[:,idx[i]], polMat, m = m)
        
        ## Construct linear operator vector
        # Compute PHS components
        for j∈1:n
            b1[j] = ϕ_xi(X, nodes[:,idx[i][j]], m, ii)
        end

        # Compute polynomial components
        for j∈1:P
            # Store polynomial term degrees and decrease derivative term
            p = polMat[:,j]
            p[ii] -= 1

            # Compute ∂/∂x of polynomial
            if p[ii] >= 0
                b2[j] = (p[ii]+1)*X[1]^p[1]*X[2]^p[2]
            else
                b2[j] = 0
            end
        end

        # Combine PHS and polynomial of linear operator vector
        b = [b1;b2]
        
        # solve for the ith node weights
        λ = A\b

        # Store weights in the discretized operator matrix
        Dx[i,idx[i]] = λ[1:n]
    end

    return Dx
end

# ∂x² discretizer 
function discretize∂xii(nodes, n, m, o, ii)
    # Number of nodes
    N = size(nodes,2)

    # Compute nearest neighbors
    idx = knnFull(nodes,n)

    # Compute polynomial degree matrix
    polMat = polyMat(o)
    
    # Number of polynomial terms
    P = size(polMat,2)

    # Preallocate space for the discretized operator matrix and the
    # linear operator vector
    Dxx = spzeros(N,N)
    b1 = zeros(n)
    b2 = zeros(P)
    
    # Begin populating discretized operator matrix
    for i∈1:N
        # Store node of interest
        X = nodes[:,i]

        # Compute collocation matrix
        A = colloc(nodes[:,idx[i]], polMat, m = m)
        
        ## Construct linear operator vector
        # Compute PHS components
        for j∈1:n
            b1[j] = ϕ_xii(X, nodes[:,idx[i][j]], m, ii)
        end

        # Compute polynomial components
        for j∈1:P
            # Store polynomial term degrees and decrease derivative term
            p = polMat[:,j]
            p[ii] -= 2

            # Compute ∂/∂x of polynomial
            if p[ii] >= 0
                b2[j] = (p[ii]+2)*(p[ii]+1)*X[1]^p[1]*X[2]^p[2]
            else
                b2[j] = 0
            end
        end

        # Combine PHS and polynomial of linear operator vector
        b = [b1;b2]
        
        # solve for the ith node weights
        λ = A\b

        # Store weights in the discretized operator matrix
        Dxx[i,idx[i]] = λ[1:n]
    end

    return Dxx
end

# ∂x∂y discretizer 
function discretize∂x∂y(nodes, n, m, o)
    # Number of nodes
    N = size(nodes,2)

    # Compute nearest neighbors
    idx = knnFull(nodes,n)

    # Compute polynomial degree matrix
    polMat = polyMat(o)
    
    # Number of polynomial terms
    P = size(polMat,2)

    # Preallocate space for the discretized operator matrix and the
    # linear operator vector
    Dxy = spzeros(N,N)
    b1 = zeros(n)
    b2 = zeros(P)
    
    # Begin populating discretized operator matrix
    for i∈1:N
        # Store node of interest
        X = nodes[:,i]

        # Compute collocation matrix
        A = colloc(nodes[:,idx[i]], polMat, m = m)
        
        ## Construct linear operator vector
        # Compute PHS components
        for j∈1:n
            b1[j] = ϕ_xij(X, nodes[:,idx[i][j]], m)
        end

        # Compute polynomial components
        for j∈1:P
            # Store polynomial term degrees and decrease derivative term
            p = polMat[:,j]
            p .-= 1

            # Compute ∂/∂x of polynomial
            if p[1] >= 0 && p[2] >= 0
                b2[j] = (p[1]+1)*(p[2]+1)*X[1]^(p[1])*X[2]^(p[2])
            else
                b2[j] = 0
            end
        end

        # Combine PHS and polynomial of linear operator vector
        b = [b1;b2]
        
        # solve for the ith node weights
        λ = A\b

        # Store weights in the discretized operator matrix
        Dxy[i,idx[i]] = λ[1:n]
    end

    return Dxy
end

## Newton's Method for finding zero level-set
function newtonSolve(X, nodes, F; n=5, m=5, o=1, ε=10^(-10), maxIts=100, prompt = false)
    # Compute number of nodes in initial guess
    N = size(X, 2)

    # Compute node bounds
    maxX = maximum(nodes[1,:])
    minX = minimum(nodes[1,:])
    maxY = maximum(nodes[2,:])
    minY = minimum(nodes[2,:])

    # Compute polynomial terms and number of terms
    polMat = polyMat(o)
    P = size(polMat,2)

    # Compute kd tree of nodes for quick nearest neighbors search
    kdtree = KDTree(nodes)

    # Initialize f*∇f/||∇f||, and function value arrays
    ∇f = zeros(2,N)
    f = zeros(N)

    # Begin Newton iteration
    for i∈1:maxIts
        # Compute nearest neighbors
        idx, tmp = knn(kdtree, X, n, true)
        
        # Compute gradient and function values at each node
        for j∈1:N
            # Compute interpolant weights
            λ = findλ(nodes[:,idx[j]], F[idx[j]], polMat, m=m)
            
            # Compute function values and first order derivatives
            f[j] = S(nodes[:,idx[j]], X[:,j], λ, polMat, m=m)
            fx = S_xi(nodes[:,idx[j]], X[:,j], λ, 1, polMat, m=m)
            fy = S_xi(nodes[:,idx[j]], X[:,j], λ, 2, polMat, m=m)

            # Compute scaled gradient
            ∇f[:,j] = f[j]*[fx,fy]/(norm([fx,fy])^2)
        end
        
        # Test convergence criteria
        if maximum(abs.(f)) < ε
            if prompt
                print(string("Surface converged in ", i, " iteration."))
            end
            return X
        end

        # Compute next solution
        X -= ∇f

        # Check bounds
        for j∈1:N
            if !(X[1,j] > minX && X[1,j] < maxX && X[2,j] > minY && X[2,j] < maxY)
                X[:,j] = nodes[:,j]
            end
        end
    end

    # Display maximum iterations used
    if prompt
        print("Maximum iterations used, solution may be inaccurate!\n")
    end
    
    return X
end
