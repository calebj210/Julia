### Parameterization Essentials

## Including used packages
using NearestNeighbors
using LinearAlgebra
using SparseArrays
using Arpack

## Backend functions
### Approximate normal functions definitions
# Inverse iteration routine to find the eigenvector associated with the
# smallest eigenvalue 
"""
    invIt(A; ϵ = 10^(-5), maxIts = 10)

Inverse iteration routine to find the eigenvector associated with the
smallest eigenvalue of the matrix `A` to a tolerance defined by `ϵ`.

`maxIts` (default = `10`) determines the max number of iterations used to find λ.

If A is singular, `invIt` returns a single vector from the nullspace of `A`.
"""
function invIt(A; ε = 10^(-5), maxIts = 10)
    # If A is singular bypass inverse iteration and return the nullspace instead
    if det(A) != 0
        B = lu(A);
        x = rand(size(A,1))
        y1 = x;
        for i ∈ 1:maxIts
            y0 = y1;
            y1 = B\x;
            x = y1/norm(y1, Inf);
            
            if norm(y0-y1, Inf) < ε || norm(y0+y1, Inf) < ε
                return x
            end
        end
    else
        return nullspace(A)[:,1]
    end

    return x
end

# Power iteration routine to find eigenvector with smallest associated
# eigenvalue.
function powIt(A; ε = 10^(-5), maxIts = 10)
    # Create random guess for the eigenvector
    x = rand(size(A,1));

    # Begin power iteration
    y1 = x;
    λ0 = 0;
    for i ∈ 1:maxIts
        λ1 = λ0
        y0 = y1;
        y1 = A*x;
        λ0 = norm(y1,Inf);
        x = y1/λ0;
        if abs(λ0-λ1) < ε
            return x
        end
    end
    
    return x
end

# Function for computing the covariance of a cluster xⱼ
"""
    covar(nodes)

Compute the covariance matrix given a cluster xⱼ defined by `nodes`
"""
function covar(nodes)
    # Number of nodes
    N = size(nodes,2);
    
    # Find average vector
    x̄ = nodes[:,1]
    for i ∈ 2:N
        x̄ += nodes[:,2];
    end
    x̄ = x̄ / N;

    # Compute covariance
    A = zeros(size(nodes,1),size(nodes,1))
    for i ∈ 1:N
        diff = nodes[:,i] - x̄;
        
        A += diff*diff';
    end
    
    return A
end


# Vector orientation algorithm
"""
    orVecs(nodes, vecs, idx)

Vector orientation algorithm that will flip all vectors defined by
`vecs` on a surface defined by `nodes` as to have the vectors all
oriented on the same side of the surface.
"""
function orVecs(nodes, vecs, idx)
    # Number of dimensions and nodes
    D,N = size(nodes);
    # Number of nodes in cluster
    n = size(idx[1]);

    # Populate index vectors
    an = repeat(1:N, inner = n);
    bn = vcat(idx...);
    D_n = spzeros(Float64, N,N);

    for i ∈ 1:size(an,1)
        ann = an[i];
        bnn = bn[i];
        D_n[ann,bnn] += vecs[:,ann] ⋅ vecs[:,bnn];
    end

    D_n = D_n'*D_n;
    # maxEig = powIt(D_n, ε = 10^(-15), maxIts = 100);
    tmp, maxEig = eigs(D_n, nev = 1, which = :LM);
    maxEig = sign.(maxEig);

    # Orient vectors
    newVecs = maxEig' .* vecs;

    return newVecs
end

# Approximate normals using associated covariance
"""
    approxNormals(nodes, idx)

Computes approximate normals given `nodes` and the nearest neighbor
indices given by `idx`.

The normals are computed using the associated covariance method. The
normal vectors are stored as column vectors in a matrix.
"""
function approxNormals(nodes, idx)
    # Number of nodes
    N = size(nodes,2);

    # Compute approximate normal by returning the eigenvector
    # associated with the smallest λ of the covariance matrix
    nmls = zeros(size(nodes,1), N);
    for i ∈ 1:N
        c = covar(nodes[:,idx[i]])
        nmls[:,i] = invIt(c, ε = 0, maxIts = 1000);
    end

    for i ∈ 1:N
        nmls[:,i] /= norm(nmls[:,i]);
    end
    
    return nmls
end
approxNormals(nodes) = approxNormals(nodes, knnFull(nodes,5))

## Functions for generating approximate normals using a sorted node
## distribution
## ((((((Depreciated normal approximator))))))
# function for finding approximate normals
# function approxNormals(nodes)
#     # Number of nodes
#     N = size(nodes,2);
#     nmls = zeros(2,N);
    
#     # First Node
#     tmp = nodes[:,2] - nodes[:,N];
#     tmp = rot90(tmp);
#     nmls[:,1] = tmp;

#     # Middle Nodes
#     for i ∈ 2:N-1
#         tmp = nodes[:,i+1] - nodes[:, i-1];
#         tmp = rot90(tmp);
#         nmls[:,i] = tmp;
#     end

#     # Last Node
#     tmp = nodes[:,1] - nodes[:,N-1];
#     tmp = rot90(tmp);
#     nmls[:,N] = tmp;

#     return nmls
# end

# Function that rotates vectors -π/2
function rot90(vec)
    tmp = zeros(2);
    tmp[1] = vec[2];
    tmp[2] = -vec[1];
    
    return tmp
end

## Find indices for k nearest neighbors
# KNN for an unordered data set
function knnFull(nodes, n)
    kdtree = KDTree(nodes);
    idx, tmp = knn(kdtree, nodes, n, true);
    
    return idx
end

# KNN for ordered data set
function knnFullOrd(nodes, n)::Array{Array{Int,1},1}
    N = size(nodes,2);
    idx = [];
    for i ∈ 1:N
        push!(idx,zeros(n))
    end
    
    # Put nth node as nth first index
    for i ∈ 1:N
        idx[i][1] = i;
    end

    # Populate rest of indices
    for i ∈ 1:N, j ∈ 2:n
        if j%2 == 0
            k = i + floor(j/2);
            idx[i][j] = k>N ? k - N : k;
        else
            k = i - floor(j/2);
            idx[i][j] = k<=0 ? N+k : k;
        end 
    end
    
    return idx
end

## enable ordered knn search
# knnFull(nodes, n) = knnFullOrd(nodes, n);

## Functions for setting up RBF envrionment
# Center main node
function center(nodes)
    nodes = nodes .- nodes[:,1];

    return nodes
end

# Finds angle between center vector and ĵ
sgn(x) = x >= 0 ? 1 : -1;
function findAngle(node)
    b = node[2]/norm(node);
    θ = sgn(node[1])acos(b);
    
    return θ
end

# Rotate nodes about origin
function rotate(nodes, θ)
    # Compute rotation matrix
    rotor = [cos(θ) -sin(θ);
             sin(θ) cos(θ)];
    
    # Rotate all vectors
    nodes = rotor*nodes;
    
    return nodes
end
