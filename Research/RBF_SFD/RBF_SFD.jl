#=
# Main functions for surface finite differences over 2D manifolds
#
# Author: Caleb Jacobs
# DLM: 2022-06-07
=#

# Needed packages and files
using LinearAlgebra
using NearestNeighbors
include("PHS.jl")

# KNN search for unordered data sets
"""
    knnFull(nodes, n)

KNN search for an unordered data set defined by `nodes` where there are `n`
neighbors
"""
function knnFull(nodes, n)
    kdtree   = KDTree(nodes)
    idx, tmp = knn(kdtree, nodes, n, true)
    
    return idx
end

# Vector orientation algorithm
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
    D_n = spzeros(Float64, N,N)

    for i ∈ 1:size(an,1)
        ann = an[i]
        bnn = bn[i]
        D_n[ann,bnn] += vecs[:,ann] ⋅ vecs[:,bnn]
    end

    D_n = D_n'*D_n
    tmp, maxEig = eigs(D_n, nev = 1, which = :LM)
    maxEig = sign.(maxEig)

    # Orient vectors
    newVecs = maxEig' .* vecs

    return newVecs
end


# Function for computing the covariance of a cluster xⱼ
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

# Approximate normals using associated covariance
"""
    appNorms(nodes, idx)

Computes approximate normals given `nodes` and the nearest neighbor
indices given by `idx`.

The normals are computed using the associated covariance method. The
normal vectors are stored as column vectors in a matrix.
"""
function appNorms(nodes, idx)
    # Number of nodes
    N = size(nodes,2)

    # Compute approximate normal by returning the eigenvector
    # associated with the smallest λ of the covariance matrix
    nmls = zeros(size(nodes,1), N)
    for i ∈ 1:N
        c = covar(nodes[:,idx[i]])
        nmls[:,i] = eigvecs(c, sortby = abs)[:, 1]
    end
    
    return nmls
end

# Node centering function given where offset is the first node in the
# array
"""
    cent(nodes)

Centers all the `nodes` so that the first node in `nodes` is at the origin.
"""
function cent(nodes)
    cents = nodes .- nodes[:,1]

    return cents
end

# Positive Nth-D rotation function based on Givens rotations
"""
    rotUp(nodes, vec)

rotate all the `nodes` so that `vec` is pointing in the posive
last coordinate.

The algorithm uses givens rotation and as such also produces a matrix of
sines and cosines so that the rotation can be reversed using `rotBack`.

`rotUp` returns a tuple of the form `(rotNodes, sc)` where sc is the
sine and cosine matrix.
"""
function rotUp(nodes, vec)
    # Dimension of space
    D = size(nodes,1)
    cluster = copy(nodes)
    nml = copy(vec)
    
    # Begin rotations
    sc = zeros(2,D-1)
    for i ∈ 1:D-1
        x   = nml[i]
        y   = nml[D]
        mag = norm([x,y])

        # Compute sine and cosine values
        s = x / mag
        c = y / mag
        
        # Store sine and cosine values for later
        sc[:,i] = [s,c]
        
        # Construct rotation matrix
        rot = [c -s;
               s c]
        
        nml[[i,D]] = [0, mag]
        cluster[[i, D], :] = rot * cluster[[i, D], :]
    end
    
    return (cluster, sc)
end

# Degree array generator
"""
    polyMat(vars, deg)

Given the number of dimensions or variables, `vars`, and the highest
degree of the polynomial specified by `deg`, polyMat produces a
polynomial power matrix.

Each column of the matrix corresponds to each term in the polynomial
and each row corresponds to the variable.

# Example
```jldoctest
julia> polyMat(2,2)
2×6 Array{Int64,2}:
 0  0  0  1  1  2
 0  1  2  0  1  0
```
This matrix will equate to the polynomial:
1 + y + y^2 + x + xy + x^2
"""
function polyMat(vars, deg)
    if deg < 0
        return Array{Int64}(undef,vars,0)
    else
        tmp::Array{Int64,1} = zeros(vars)
        array::Array{Array{Int64,1},1} = [copy(tmp)]
        idx = 1
        i = 1
        while true
            if tmp[1] == deg
                break
            end
            if sum(tmp) == deg
                tmp[idx] = 0
                tmp[idx-1] += 1
                idx -= 1
                push!(array,copy(tmp))
            elseif idx != vars
                idx += 1
            else
                tmp[idx] += 1
                push!(array,copy(tmp))
            end
        end
        return reduce(hcat, array)
    end
end

