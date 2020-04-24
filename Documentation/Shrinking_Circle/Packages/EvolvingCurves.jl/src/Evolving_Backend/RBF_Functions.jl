## RBF related functions

# RBF definitons
# 1D
# Gaussian
# ϕ(x1,x2,m) = exp(-abs2(m*(x1-x2)));
# ϕ_x(x1,x2,m) = 2*(m^2)*(x2-x1)*exp(-abs2(m*(x1-x2)));

# PHS
ϕ(x1,x2,m) = abs(x1-x2)^m;
ϕ_x(x1,x2,m) = m*sign(x1-x2)*abs(x1-x2)^(m-1);
ϕ_xx(x1,x2,m) = m*(m-1)*abs(x1-x2)^(m-2);
ϕ_xxx(x1,x2,m) = m*(m-1)*(m-2)*sign(x1-x2)*abs(x1-x2)^(m-3);
ϕ_xxxx(x1,x2,m) = m*(m-1)*(m-2)*(m-3)*abs(x1-x2)^(m-4);

# Collocation matrix builder
function collocM(x, m, o)
    n = size(x, 1);
    A0 = zeros(n,n);
    A1 = zeros(o+1,n);

    # Define A without helping terms
    for i ∈ 1:n, j ∈ 1:n
        A0[i,j] = ϕ(x[i], x[j], m);
    end
    
    # Create helping terms
    A1[1,:] .= 1;
    for i ∈ 1:o, j ∈ 1:n
        A1[i+1,j] = x[j]^i;
    end

    # Create full A
    A = [A0 A1';
         A1 zeros(o+1,o+1)];

    return A
end

# Interpolation weights for 1D with helping terms
function interpolate(nodes, m, o)
    n = size(nodes, 2);

    # Compute collocation matrix with polynomial terms
    A = collocM(nodes[1,:],m,o);
    
    # Store local function values
    f = [nodes[2,:]; zeros(o+1)];
    
    # Solve for weights
    λ = A\f;

   return λ
end

# Interpolation function
function S(t, x, λ, m)
    n = size(x,1);
    nn = size(λ,1) - n;
    nt = size(t,1);
    
    s = zeros(nt);
    for i ∈ 1:n, j ∈ 1:nt
        s[j] += λ[i]*ϕ(t[j], x[i], m);
    end
    for i ∈ 0:nn-1, j ∈ 1:nt
        s[j] += λ[n+i+1]*t[j]^i;
    end
    
    return s
end 

# RBF interpolant derivative at s=0
function S_x(x, λ, m)
    n = size(x,1);
    
    s = 0;
    for i ∈ 1:n
        s += λ[i]*ϕ_x(0, x[i], m);
    end

    # Add polynomial contribution if there is one
    if size(λ,1) > n+1
        s += λ[n+2];
    end
    
    return s
end

# RBF interpolant second derivative at s=0
function S_xx(x, λ, m)
    n = size(x,1);
    
    s = 0;
    for i ∈ 1:n
        s += λ[i]*ϕ_xx(0, x[i], m);
    end

    # Add polynomial contribution if there is one
    if size(λ,1) > n+2
        s += 2*λ[n+3];
    end
    
    return s
end

# RBF interpolant third derivative at s=0
function S_xxx(x, λ, m)
    n = size(x,1);
    
    s = 0;
    for i ∈ 1:n
        s += λ[i]*ϕ_xxx(0, x[i], m);
    end

    # Add polynomial contribution if there is one
    if size(λ,1) > n+3
        s += 6*λ[n+4];
    end
    
    return s
end

# RBF interpolant third derivative at s=0
function S_xxxx(x, λ, m)
    n = size(x,1);
    
    s = 0;
    for i ∈ 1:n
        s += λ[i]*ϕ_xxxx(0, x[i], m);
    end

    # Add polynomial contribution if there is one
    if size(λ,1) > n+4
        s += 24*λ[n+5];
    end
    
    return s
end
