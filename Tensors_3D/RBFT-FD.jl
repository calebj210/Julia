### RBFT-FD functions

## Including used packages
using NearestNeighbors
using LinearAlgebra
using SparseArrays
using Arpack
include("PHS.jl")

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
function invIt(A; ϵ = 10^(-5), maxIts = 10)
    # If A is singular bypass inverse iteration and return the nullspace instead
    if det(A) != 0
        B = lu(A);
        x = rand(size(A,1))
        y1 = x;
        for i ∈ 1:maxIts
            y0 = y1;
            y1 = B\x;
            x = y1/norm(y1, Inf);
            
            if norm(y0-y1, Inf) < ϵ || norm(y0+y1, Inf) < ϵ
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

# KNN search for unordered data sets
"""
    knnFull(nodes, n)

KNN search for an unordered data set defined by `nodes` where there are `n`
neighbors
"""
function knnFull(nodes, n)
    kdtree = KDTree(nodes);
    idx, tmp = knn(kdtree, nodes, n, true);
    
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
    tmp, maxEig = eigs(D_n, nev = 1, which = :LM);
    maxEig = sign.(maxEig);

    # Orient vectors
    newVecs = maxEig' .* vecs;

    return newVecs
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
    N = size(nodes,2);

    # Compute approximate normal by returning the eigenvector
    # associated with the smallest λ of the covariance matrix
    nmls = zeros(size(nodes,1), N);
    for i ∈ 1:N
        c = covar(nodes[:,idx[i]])
        nmls[:,i] = invIt(c, ϵ = 0, maxIts = 1000);
    end
    
    return nmls
end


## Node orienting functions
# Node centering function given where offset is the first node in the
# array
"""
    cent(nodes)

Centers all the `nodes` so that the first node in `nodes` is at the origin.
"""
function cent(nodes)
    cents = nodes .- nodes[:,1];

    return cents
end

# Alternative center function that replaces given nodes
""" 
    cent!(nodes)

same as `cent` but stores the centered nodes in `nodes`.
"""
function cent!(nodes)
    nodes .-= nodes[:,1];
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
    D= size(nodes,1);
    cluster = copy(nodes);
    nml = copy(vec);
    
    # Begin rotations
    sc = zeros(2,D-1)
    for i ∈ 1:D-1
        x = nml[i];
        y = nml[D];
        mag = norm([x,y])

        # Compute sine and cosine values
        s = x/mag;
        c = y/mag;
        
        # Store sine and cosine values for later
        sc[:,i] = [s,c];
        
        # Construct rotation matrix
        rot = [c -s;
               s c];
        
        nml[[i,D]] = [0,mag]
        cluster[[i,D],:] = rot*cluster[[i,D],:];
    end
    
    return (cluster,sc)
end

# Alternative rotUp function that replaces given node sets
"""
    rotUp!(nodes, normal)

Same as `rotUp` but now stores the rotated nodes in `nodes` and only
returns the sine and cosine matrix.
"""
function rotUp!(nodes, normal)
    # Dimension of space
    D= size(nodes,1);
    nml = copy(normal);
    
    # Begin rotations
    sc = zeros(2,D-1)
    for i ∈ 1:D-1
        x = nml[i];
        y = nml[D];
        mag = norm([x,y])

        # Compute sine and cosine values
        s = x/mag;
        c = y/mag;
        
        # Store sine and cosine values for later
        sc[:,i] = [s,c];
        
        # Construct rotation matrix
        rot = [c -s;
               s c];
        
        nml[[i,D]] = [0,mag]
        nodes[[i,D],:] = rot*nodes[[i,D],:];
    end
    
    return sc
end

# Reverse rotation function given sine and cosines of previous rotations
"""
    rotBack(vector, sc)

Undo the rotations of `rotUp` on single vector given a sine and cosine
matrix, `sc`
"""
function rotBack(vector, sc)
    D = size(vector,1);
    vec = copy(vector);
    
    for i ∈ D-1:-1:1
        s = sc[1,i];
        c = sc[2,i];
        
        # Construct rotation matrix
        rot = [c s;
               -s c];
        
        vec[[i,end]] = rot*vec[[i,end]];
    end

    return vec
end

# Alternative rotBack function that replaces vector
"""
    rotBack!(vector, sc)

Same as `rotBack` but now stores the rotated vector in `vector`
"""
function rotBack!(vector, sc)
    D = size(vector,1);
    
    for i ∈ D-1:-1:1
        s = sc[1,i];
        c = sc[2,i];
        
        # Construct rotation matrix
        rot = [c s;
               -s c];
        
        vector[[i,end]] = rot*vector[[i,end]];
    end
end

## Polynomial functions
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
        tmp::Array{Int64,1} = zeros(vars);
        array::Array{Array{Int64,1},1} = [copy(tmp)];
        idx = 1;
        i = 1;
        while true
            if tmp[1] == deg
                break
            end
            if sum(tmp) == deg
                tmp[idx] = 0;
                tmp[idx-1] += 1;
                idx -= 1;
                push!(array,copy(tmp));
            elseif idx != vars
                idx += 1;
            else
                tmp[idx] += 1;
                push!(array,copy(tmp));
            end
        end
        return reduce(hcat, array)
    end
end

## Finite difference functions
# Collocation matrix generator
"""
    colloc(nodes, poly = Array{Int64}(undef,0,0); m = 3)

`colloc` takes the local nodes, `nodes`, and polynomial terms, `poly`, to
produce a collocation matrix for RBF interpolation or RBF-FD operator
discretization. The order of the PHS can specified by `m`.
"""
function colloc(nodes, poly = Array{Int64}(undef,0,0); m = 3)
    # Number of nodes in cluster
    N = size(nodes,2);
    # Number of dimensions
    D = size(nodes,1);
    # Number of polynomial terms
    P = size(poly,2);
    
    # Preallocating space
    A11 = zeros(N,N);
    A12 = zeros(N,P);

    # Add PHS contributions
    for j ∈ 1:N, i ∈ 1:N
        A11[i,j] = ϕ(nodes[:,j],nodes[:,i],m);
    end

    # Add polynomial contributions
    for j ∈ 1:P, i ∈ 1:N
        tmp = 1;
        for k ∈ 1:D
            tmp *= nodes[k,i]^poly[k,j];
        end
            
        A12[i,j] = tmp;
    end

    # Construst A
    A = [A11 A12;
         A12' zeros(P,P)];
    
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
function findλ(nodes, f, poly = Array{Int64}(undef,0,0); m = 3)
    # Number of nodes
    N = size(f,1);
    # Number of dimensions
    D = size(nodes,1);
    # Number of polynomial terms
    P = size(poly,2);

    # Construct collocation matrix
    A = colloc(nodes, poly, m = m);

    # Construct b in Aλ=b
    b = [f;zeros(P)];

    # Solve for the weights
    λ = A\b;

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
    N = size(nodes,2);
    # Number of dimensions
    D = size(nodes,1);
    # Number of polynomial terms
    P = size(poly,2);

    # Add PHS contribution
    tmp = 0;
    for i ∈ 1:N
        tmp += λ[i]*ϕ(Xc,nodes[:,i],m);
    end

    # Add polynomial contribution
    for i ∈ 1:P
        tmp2 = λ[i+N];
        for j ∈ 1:D
            tmp2 *= Xc[j]^poly[j,i];
        end
        tmp += tmp2
    end

    return tmp
end

# Compute derivative of interpolant at X_c with respect to xi
"""
    S_xi(nodes, Xc, λ, ii, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)

Compute the derivative of the RBF interpolation at `Xc` with respect to the
`ii`th coordinate given the weights `λ`.

See `S` or `findλ`.
"""
function S_xi(nodes, Xc, λ, ii, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)
    # Number of nodes
    N = size(nodes,2);
    # Number of dimensions
    D = size(nodes,1);
    # Number of polynomial terms
    P = size(poly,2);

    # Add PHS contribution
    tmp = 0;
    for i ∈ 1:N
        tmp += λ[i]*ϕ_xi(Xc,nodes[:,i],m,ii);
    end
    
    # Add polynomial contribution
    for i ∈ 1:P
        p = poly[ii,i];
        tmp2 = 0;
        if p-1 >= 0
            tmp2 = λ[i+N]*p;
            for j ∈ 1:ii-1
                tmp2 *= Xc[j]^poly[j,i];
            end

            # Add derivative component of polynomial
            tmp2 *= Xc[ii]^(p-1);
            
            for j ∈ ii+1:D
                tmp2 *= Xc[j]^poly[j,i];
            end
        end
        
        tmp += tmp2
    end

    return tmp
end

# Compute derivative of interpolant at X_c with respect to xi
"""
    S_xii(nodes, Xc, λ, ii, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)

Compute the second derivative of the RBF interpolation at `Xc` with respect to the
`ii`th coordinate given the weights `λ`.

See `S`, `S_xi`, or `findλ`.
"""
function S_xii(nodes, Xc, λ, ii, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)
    # Number of nodes
    N = size(nodes,2);
    # Number of dimensions
    D = size(nodes,1);
    # Number of polynomial terms
    P = size(poly,2);

    # Add PHS contribution
    tmp = 0;
    for i ∈ 1:N
        tmp += λ[i]*ϕ_xii(Xc,nodes[:,i],m,ii);
    end
    
    # Add polynomial contribution
    for i ∈ 1:P
        p = poly[ii,i];
        tmp2 = 0;
        if p-2 >= 0
            tmp2 = λ[i+N]*p*(p-1);
            
            # Add standard polynomial terms
            for j ∈ 1:ii-1
                tmp2 *= Xc[j]^poly[j,i]
            end
            
            # Add derivative component of polynomial
            p = p-2;
            tmp2 *= Xc[ii]^p;
            
            # Add more standard polynomial terms
            for j ∈ ii+1:D
                tmp2 *= Xc[j]^poly[j,i]
            end
        end
        
        tmp += tmp2
    end

    return tmp
end

# Compute derivative of interpolant at X_c with respect to xi
"""
    S_xii(nodes, Xc, λ, ii, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)

Compute the mixed partial derivative of the RBF interpolation at `Xc` with respect to the
`ii`th and `jj`th coordinates given the weights `λ`.

See `S`, `S_xi`, `S_xij`, or `findλ`.
"""
function S_xij(nodes, Xc, λ, ii, jj, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)
    # Number of nodes
    N = size(nodes,2);
    # Number of dimensions
    D = size(nodes,1);
    # Number of polynomial terms
    P = size(poly,2);

    # Add PHS contribution
    tmp = 0;
    for i ∈ 1:N
        tmp += λ[i]*ϕ_xii(Xc,nodes[:,i],m,ii);
    end
    
    # Add polynomial contribution
    for i ∈ 1:P
        (k1,k2) = (ii<jj ? (ii,jj) : (jj,ii));
        p1 = poly[k1,i];
        p2 = poly[k2,i];
        tmp2 = λ[i+N]*p1*p2;
        
        # Add standard polynomial terms
        for j ∈ 1:k1-1
            tmp2 *= Xc[j]^poly[j,i]
        end
        
        # Add derivative component of polynomial
        p1 = p1-1;
        if p1 >= 0
            tmp2 *= Xc[ii]^p1;
        end
        
        for j ∈ k1+1:k2-1
            tmp2 *= Xc[j]^poly[j,i]
        end
        
        # Add derivative component of polynomial
        p2 = p2-1;
        if p2 >= 0
            tmp2 *= Xc[ii]^p2;
        end

        # Add more standard polynomial terms
        for j ∈ k2+1:D
            tmp2 *= Xc[j]^poly[j,i]
        end
        
        tmp += tmp2
    end

    return tmp
end

# Function for discretizing the Laplace-Beltrami operator on a 2-D
# surface embedded in 3-D ambient space
function lapDisc3D(nodes, n, deg, m)
    # Number of dimensions and nodes
    (D,N) = size(nodes);

    idx = knnFull(nodes,n);

    appnms = appNorms(nodes, idx);
end
