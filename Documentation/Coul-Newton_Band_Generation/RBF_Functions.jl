## RBF related functions

# PHS
ϕ(X,Y,m) = norm(X-Y)^m;
ϕ_xi(X,Y,m,i) = m*(X[i]-Y[i])*norm(X-Y)^(m-2)
ϕ_xii(X,Y,m,i) =
    begin
        nm = norm(X-Y);
        
        return m*nm^(m-2) + m*(m-2)*(X[i]-Y[i])^2*nm^(m-4)
    end
ϕ_xij(X,Y,m) = m*(m-2)*(X[1]-Y[1])*(X[2]-Y[2])*norm(X-Y)^(m-4);

# Polynomial degree array generator
"""
    polyMat(deg)

Given highest degree of the polynomial specified by `deg`, polyMat
produces a polynomial power matrix.

Each column of the matrix corresponds to each term in the polynomial
and each row corresponds to the variable.

# Example
```jldoctest
julia> polyMat(2)
2×6 Array{Int64,2}:
 0  0  0  1  1  2
 0  1  2  0  1  0
```
This matrix will equate to the polynomial:
1 + y + y^2 + x + xy + x^2
"""
function polyMat(deg)
    if deg < 0
        return Array{Int64}(undef,2,0)
    else
        tmp = zeros(2);
        array = [copy(tmp)];
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
            elseif idx != 2
                idx += 1;
            else
                tmp[idx] += 1;
                push!(array,copy(tmp));
            end
        end
        return reduce(hcat, array)
    end
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
    N = size(nodes,2);
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
        for k ∈ 1:2
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
        for j ∈ 1:2
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
            
            for j ∈ ii+1:2
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
        p = poly[ii,i]
        tmp2 = 0
        if p-2 >= 0
            tmp2 = λ[i+N]*p*(p-1)
            
            # Add standard polynomial terms
            for j ∈ 1:ii-1
                tmp2 *= Xc[j]^poly[j,i]
            end
            
            # Add derivative component of polynomial
            tmp2 *= Xc[ii]^(p-2)
            
            # Add more standard polynomial terms
            for j ∈ ii+1:2
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
function S_xij(nodes, Xc, λ, poly = Array{Int64}(undef,size(nodes,1),0); m = 3)
    # Number of nodes
    N = size(nodes,2);
    # Number of dimensions
    D = size(nodes,1);
    # Number of polynomial terms
    P = size(poly,2);

    # Add PHS contribution
    tmp = 0;
    for i ∈ 1:N
        tmp += λ[i]*ϕ_xij(Xc,nodes[:,i],m);
    end
    
    # Add polynomial contribution
    for i ∈ 1:P
        polip = poly[:,i]

        tmp2 = 0.0
        if polip[1]-1 >= 0 && polip[2]-1 >= 0
            tmp2 = λ[N+i]*polip[1]*polip[2]
            polip .-= 1

            for j = 1:2
                tmp2 *= Xc[j]^polip[j]
            end
        end
        tmp += tmp2
    end

    return tmp
end
