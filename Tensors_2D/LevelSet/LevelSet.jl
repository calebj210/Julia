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
function newtonSolve(X, nodes, F; n=5, m=5, o=1, ε=10^(-10), maxIts=100)
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
            return X
        end

        # Compute next solution
        X -= .9*∇f

        # Check bounds
        for j∈1:N
            if !(X[1,j] > minX && X[1,j] < maxX && X[2,j] > minY && X[2,j] < maxY)
                X[:,j] = nodes[:,j]
            end
        end
    end

    # Display maximum iterations used
    print("Maximum iterations used, solution may be inaccurate!\n")
    return X
end
