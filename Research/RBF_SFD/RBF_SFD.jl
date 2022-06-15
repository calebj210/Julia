#=
# Main functions for surface finite differences over 2D manifolds
#
# Author: Caleb Jacobs
# DLM: 2022-06-10
=#

# Needed packages and files
using LinearAlgebra
using NearestNeighbors
include("PHS.jl")
include("SurfaceElements.jl")

# Handy data type for RBF-SFD
struct Common
    λs
    idx
    pols
    clusts
    scs
    n
    m
    o
end


# KNN search for unordered data sets
"""
    knnFull(nodes, n)

KNN search for an unordered data set defined by `nodes` where there are `n`
neighbors
"""
function knnFull(nodes::Matrix{Float64}, n)
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

rotate all the `nodes` so that `vec` is pointing in the positive
last coordinate.
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
        mag = hypot(x, y)
        
        if mag == 0
            sc[:, i] = [0, 0]
            continue
        end

        # Compute sine and cosine values
        s = x / mag
        c = y / mag

        sc[:, i] = [s, c]
        
        # Construct rotation matrix
        rot = [c -s;
               s  c]
        
        nml[[i,D]] = [0, mag]
        cluster[[i, D], :] = rot * cluster[[i, D], :]
    end
    
    return (cluster, sc)
end

# Reverse rotation function given sine and cosines of previous rotations
"""
    rotBack(vector, sc)

Undo the rotations of `rotUp` on single vector given a sine and cosine
matrix, `sc`
"""
function rotBack(vector, sc)
    D = size(sc, 2);
    vec = copy(vector);
    
    for i ∈ D:-1:1
        s = sc[1,i];
        c = sc[2,i];
        
        # Construct rotation matrix
        rot = [c s;
               -s c];
        
        vec[[i,end]] = rot*vec[[i,end]];
    end

    return vec
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

# Compute single term in polynomial
function polyVal(x, pols)
    reduce(*, x .^ pols)
end

# Collocation matrix generator
"""
    colloc(nodes, poly = Array{Int64}(undef,0,0); m = 3)

`colloc` takes the local nodes, `nodes`, and polynomial terms, `poly`, to
produce a collocation matrix for RBF interpolation or RBF-FD operator
discretization. The order of the PHS can specified by `m`.
"""
function colloc(nodes, poly = Array{Int64}(undef,0,0); m = 3)
    # Number of nodes in cluster
    N = size(nodes,2)
    # Number of dimensions
    D = size(nodes,1)
    # Number of polynomial terms
    P = size(poly,2)
    
    # Preallocating space
    A11 = zeros(N,N)
    A12 = zeros(N,P)

    # Add PHS contributions
    for j ∈ 1:N, i ∈ 1:N
        A11[i,j] = ϕ(nodes[:,j],nodes[:,i],m)
    end

    # Add polynomial contributions
    for j ∈ 1:P, i ∈ 1:N
        tmp = 1
        for k ∈ 1:D
            tmp *= nodes[k,i]^poly[k,j]
        end
            
        A12[i,j] = tmp
    end

    # Construst A
    A = [A11 A12;
         A12' zeros(P,P)]
    
    return A 
end

# Routine for computing RBF-PHS interpolant weights
"""
    findλ(nodes, f, poly = Array{Int64}(undef,0,0); m = 3)

`findλ` computes the RBF-PHS + polynomial term interpolation weights
given the local coordinate `nodes` and the function values, `f`. `findλ`
also if no polynomial power matrix, `poly`, is given, no polynomial
terms are used. The order of the PHS, `m`, can also be specified.
"""
function findλ(nodes, poly = Array{Int64}(undef,0,0); m = 3)
    # Number of nodes
    d, N = size(nodes)

    # Number of polynomial terms
    P = size(poly,2)

    # Construct collocation matrix
    A = colloc(nodes[1:d - 1, :], poly, m = m)

    # Construct b in Aλ=b
    b = [nodes[d, :]; zeros(P)]

    # Solve for the weights
    λ = A \ b

    return λ
end

# Compute interpolant at X_c with respect to xi
"""
    S(nodes, Xc, λ, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)

Compute the RBF interpolation at `Xc` given the interpolation weights `λ`.

`λ` can be found using `findλ`.
"""
function S(nodes, Xc, λ, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)
    # Number of nodes
    N = size(nodes,2)
    # Number of dimensions
    D = size(nodes,1) - 1
    # Number of polynomial terms
    P = size(poly,2)

    # Add PHS contribution
    tmp = 0
    for i ∈ 1:N
        tmp += λ[i] * ϕ(Xc, nodes[1:D, i], m)
    end

    # Add polynomial contribution
    for i ∈ 1:P
        tmp2 = λ[i+N]
        for j ∈ 1:D
            tmp2 *= Xc[j]^poly[j,i]
        end
        tmp += tmp2
    end

    return tmp
end

# Compute commonly needed objects
function getCommons(nodes, n, m, o)
    d, N = size(nodes)                          # Dimension and # of nodes
    
    pols = polyMat(d - 1, o)                    # Polynomial matrix

    idx = knnFull(nodes, n)                     # Nearest neighbor indexing

    aNorms = appNorms(nodes, idx)               # Approximate normals
    aNorms = orVecs(nodes, aNorms, idx)         # Orient normals

    λs     = Array{Array{Float64, 1}, 1}(undef, N)
    clusts = Array{Array{Float64, 2}, 1}(undef, N)
    scs    = Array{Array{Float64, 2}, 1}(undef, N)
    for i ∈ 1:N
        clust = cent(nodes[:, idx[i]])          # Center node cluster
        clust, sc = rotUp(clust, aNorms[:, i])  # Rotate node cluster

        λs[i]     = findλ(clust, pols, m = m)   # Compute interpolant weights
        clusts[i] = clust
        scs[i]    = sc
    end 

    return Common(λs, idx, pols, clusts, scs, n, m, o)
end

# Compute normals to each node
function getNormals(nodes, c::Common)
    d, N = size(nodes)                          # Dimension and # of nodes
    
    clusts = c.clusts                           # Rotated node clusters
    pols   = c.pols                             # Polynomial matrix
    λs     = c.λs                               # Interpolant weights
    scs    = c.scs                              # rotation matrices
    m      = c.m                                # Shape parameter

    normals = zeros(d, N)                       # Allocate space for normals
    for i ∈ 1:N
        nml = normal(zeros(d - 1), x -> 
                     S(clusts[i], x, λs[i], pols, m = m))
        normals[:, i] = rotBack(nml, scs[i])
    end

    return normals
end 

# Compute mean curvature
function getH(nodes, c::Common)
    d, N = size(nodes)                          # Dimension and # of nodes
    
    λs     = c.λs                               # Interpolant weights
    pols   = c.pols                             # Polynomial matrix
    clusts = c.clusts                           # Rotated node clusters
    m      = c.m                                # Shape parameter

    Hs = zeros(N)
    for i ∈ 1:N
        Hs[i] = H(zeros(d - 1), x -> 
                     S(clusts[i], x, λs[i], pols, m = m))
    end

    return Hs
end

#### Discretize Laplace-Beltrami operator BROKEN...
function discΔ(nodes, c::Common)
    d, N = size(nodes)                          # Dimension and # of nodes
    clusts = c.clusts                           # Rotated node clusters
    pols   = c.pols                             # Polynomial matrix
    λs     = c.λs                               # Interpolant weights
    n      = c.n                                # # of nearest neighbors
    m      = c.m                                # Shape parameter
    idx    = c.idx                              # KNN indices
    P      = size(pols, 2)                      # # of polynomial terms
    
    DΔ = zeros(N, N)                            # Initialize differentiation matrix
    for i ∈ 1:N
        s0 = zeros(d - 1)                       # Center to evaluate at
        
        A = colloc(clusts[i][1:d - 1, :],       # Collocation matrix
                   pols, m = m)

        Lϕ = zeros(n + P)                       # RHS 1
        for j ∈ 1:n
            Lϕ[j] = Δ(s0,
                      x -> S(clusts[i], x, λs[i], pols, m = m),
                      x -> ϕ(x, clusts[i][1:d - 1, j], m))
        end

        for j ∈ 1:P
            Lϕ[n + j] = Δ(s0,
                      x -> S(clusts[i], x, λs[i], pols, m = m), 
                      x -> polyVal(x, pols[j]))
        end

        w = A \ Lϕ                              # Discretization weights

        DΔ[i, idx[i]] = w[1:n]                  # Build discretization matrix
    end

    return DΔ
end


# Compute Laplace-Beltrami operator
function computeΔ(nodes, F, c::Common)
    d, N = size(nodes)                          # Dimension and # of nodes
    clusts = c.clusts                           # Rotated node clusters
    pols   = c.pols                             # Polynomial matrix
    λs     = c.λs                               # Interpolant weights
    n      = c.n                                # # of nearest neighbors
    m      = c.m                                # Shape parameter
    idx    = c.idx                              # KNN indices
    P      = size(pols, 2)                      # # of polynomial terms
    
    ΔF = zeros(N)
    for i ∈ 1:N
        s0 = zeros(d - 1)                       # Center to evaluate at
        
        A = colloc(clusts[i][1:d - 1, :],       # Collocation matrix
                   pols, m = m)

        λF = A \ [F[idx[i]]; zeros(P)]

        ΔF[i] = Δ(s0,
                  x -> S(clusts[i], x, λs[i], pols, m = m),
                  x -> S(clusts[i], x, λF, pols, m = m))
    end

    return ΔF
end

# Compute hyperviscoisty term
function computeHyperV(nodes, F, c::Common)
    d, N = size(nodes)                          # Dimension and # of nodes
    clusts = c.clusts                           # Rotated node clusters
    pols   = c.pols                             # Polynomial matrix
    λs     = c.λs                               # Interpolant weights
    n      = c.n                                # # of nearest neighbors
    m      = c.m                                # Shape parameter
    idx    = c.idx                              # KNN indices
    P      = size(pols, 2)                      # # of polynomial terms
    
    Δ⁸F = zeros(N)
    for i ∈ 1:N
        s0 = zeros(d - 1)                       # Center to evaluate at
        
        A = colloc(clusts[i][1:d - 1, :],       # Collocation matrix
                   pols, m = m)

        λF = A \ [F[idx[i]]; zeros(P)]

        Δ⁸F[i] = Δ(s0,
                  x -> S(clusts[i], x, λs[i], pols, m = m),
                  x -> S(clusts[i], x, λF, pols, m = m))
    end

    return Δ⁸F
end

# Compute surface divergence term
function computeSurfDiv(nodes, F, c::Common, ε, δ)
    d, N = size(nodes)                          # Dimension and # of nodes
    clusts = c.clusts                           # Rotated node clusters
    pols   = c.pols                             # Polynomial matrix
    λs     = c.λs                               # Interpolant weights
    n      = c.n                                # # of nearest neighbors
    m      = c.m                                # Shape parameter
    idx    = c.idx                              # KNN indices
    P      = size(pols, 2)                      # # of polynomial terms
    
    Div = zeros(N)
    for i ∈ 1:N
        s0 = zeros(d - 1)                       # Center to evaluate at
        
        A = colloc(clusts[i][1:d - 1, :],       # Collocation matrix
                   pols, m = m)

        λF = A \ [F[idx[i]]; zeros(P)]

        u = x -> S(clusts[i], x, λF,    pols, m = m)        
        z = x -> S(clusts[i], x, λs[i], pols, m = m)

        Div[i] = ∇Γdot(s0,
                  z,
                  x -> (ε * H(x, z) + δ * u(x)) * normal(x, z))
    end

    return Div
end

# Compute surface divergence of u
function computeSurfDivU(nodes, F, c::Common)
    d, N = size(nodes)                          # Dimension and # of nodes
    clusts = c.clusts                           # Rotated node clusters
    pols   = c.pols                             # Polynomial matrix
    λs     = c.λs                               # Interpolant weights
    n      = c.n                                # # of nearest neighbors
    m      = c.m                                # Shape parameter
    idx    = c.idx                              # KNN indices
    P      = size(pols, 2)                      # # of polynomial terms
    
    Div = zeros(N)
    for i ∈ 1:N
        s0 = zeros(d - 1)                       # Center to evaluate at
        
        A = colloc(clusts[i][1:d - 1, :],       # Collocation matrix
                   pols, m = m)

        λF = A \ [F[idx[i]]; zeros(P)]

        u = x -> S(clusts[i], x, λF,    pols, m = m)        
        z = x -> S(clusts[i], x, λs[i], pols, m = m)

        Div[i] = ∇Γdot(s0, z, x -> u(x) * normal(x, z))
    end

    return Div
end
# Compute normal velocity
function computeNormVel(nodes, u, c::Common, ε, δ)
    d, N = size(nodes)                          # Dimension and # of nodes
    clusts = c.clusts                           # Rotated node clusters
    pols   = c.pols                             # Polynomial matrix
    λs     = c.λs                               # Interpolant weights
    n      = c.n                                # # of nearest neighbors
    m      = c.m                                # Shape parameter
    idx    = c.idx                              # KNN indices
    P      = size(pols, 2)                      # # of polynomial terms
    
    V = zeros(N)
    for i ∈ 1:N
        s0 = zeros(d - 1)                       # Center to evaluate at
        
        A = colloc(clusts[i][1:d - 1, :],       # Collocation matrix
                   pols, m = m)

        λF = A \ [u[idx[i]]; zeros(P)]
        
        z = x -> S(clusts[i], x, λs[i], pols, m = m)

        V[i] = ε * H(s0, z) + δ * u[i]
    end

    return V
end
