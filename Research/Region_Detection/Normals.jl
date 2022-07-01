#=
# Useful algorithms for computing and orienting normals to a surface
#
# Author: Caleb Jacobs
# DLM: 27-06-2022
=#

"""
    orVecs(nodes, vecs, idx)

Vector orientation algorithm that will flip all vectors defined by
`vecs` on a surface defined by `nodes` as to have the vectors all
oriented on the same side of the surface.
"""
function orVecs(nodes, vecs, idx)
    D,N = size(nodes)   # Dimensions and number of nodes
    n = size(idx[1])    # Number of nodes in cluster

    # Populate index vectors
    an = repeat(1:N, inner = n)
    bn = vcat(idx...)
    D_n = zeros(N, N)

    for i ∈ 1:size(an,1)
        ann = an[i]
        bnn = bn[i]
        D_n[ann,bnn] += vecs[:,ann] ⋅ vecs[:,bnn]
    end

    D_n = D_n'*D_n
    maxEig = eigvecs(D_n, sortby = x -> -abs(x))[:, 1]
    maxEig = sign.(maxEig)

    # Orient vectors
    newVecs = maxEig' .* vecs

    testVec = 100 * ones(D)
    outerIdx, tmp = nn(KDTree(nodes), testVec)
    testVec -= nodes[:, outerIdx[1]]

    if newVecs[:, outerIdx[1]] ⋅ testVec < 0
        newVecs *= -1
    end

    return newVecs
end

"""
    covar(nodes)

Compute the covariance matrix given a cluster xⱼ defined by `nodes`
"""
function covar(nodes)
    # Number of nodes
    N = size(nodes, 2)
    
    # Find average vector
    x̄ = nodes[:, 1]
    for i ∈ 2 : N
        x̄ += nodes[:, 2]
    end
    x̄ = x̄ / N

    # Compute covariance
    A = zeros(size(nodes, 1), size(nodes, 1))
    for i ∈ 1 : N
        diff = nodes[:, i] - x̄
        
        A += diff * diff'
    end
    
    return A
end

"""
    getn⃗(nodes, idx)

Computes approximate normals given `nodes` and the nearest neighbor
indices given by `idx`.

The normals are computed using the associated covariance method. The
normal vectors are stored as column vectors in a matrix.
"""
function getn⃗(nodes, idx)
    # Number of nodes
    N = size(nodes,2)

    # Compute approximate normal by returning the eigenvector
    # associated with the smallest λ of the covariance matrix
    nmls = zeros(size(nodes,1), N)
    for i ∈ 1:N
        c = covar(nodes[:,idx[i]])
        nmls[:,i] = eigvecs(c, sortby = abs)[:, 1]
    end

    nmls = orVecs(nodes, nmls, idx)     # Orient normals
    
    return nmls
end
