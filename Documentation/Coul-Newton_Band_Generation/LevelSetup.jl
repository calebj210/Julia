### Essential Level Set Functions

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

## Find indices for k nearest neighbors
# KNN for an unordered data set
function knnFull(nodes, n)
    kdtree = KDTree(nodes);
    idx, tmp = knn(kdtree, nodes, n, true);
    
    return idx
end
