### Steady State 1D Dirichlet B.C. Problem Solver Using
### Spectral-Chebyshev Method

using LinearAlgebra
using Plots
gr()

# Function for populating boundary rows of the collocation matrix
function rowPop(x,M)
    T = zeros(M+2);

    T[1] = 1;
    T[2] = x;
    for i ∈ 3:M+2;
        T[i] = 2x*T[i-1] - T[i-2];
    end

    return T
end

# Function for populating DE rows of the collocation matrix
function rowPop_xx(x,M)
    T_x = zeros(M+1);

    T_x[1] = 0
    T_x[2] = 1;
    T_x[3] = 2x;
    for i ∈ 4:M+1
        T_x[i] = 2x*T_x[i-1] - T_x[i-2];
    end
    T_x .*= 0:M;

    T_xx = zeros(M+2);
    for i ∈ 3:M+2
        T_xx[i] = 4T_x[i-1] + 2x*T_xx[i-1] - T_xx[i-2];
    end

    return T_xx
end

# Collocation populator
function collocA(x)
    M = size(x,1) - 2;

    A = zeros(M+2,M+2);
    A[1,:] = rowPop(x[1],M);
    A[M+2,:] = rowPop(x[M+2],M);

    for i ∈ 2:M+1
        A[i,:] = rowPop_xx(x[i],M);
    end

    return A
end

# F populator
function popF(x,α,β,f)
    M = size(x,1) - 2;

    F = zeros(M+2);
    F[1] = α;
    F[M+2] = β;
    for i ∈ 2:M+1
        F[i] = f(x[i])
    end

    return F
end

# ODE Solver
function uSolve(a,b,α,β,f,M)
    x = zeros(M+2);
    for i ∈ 0:M+1
        x[i+1] = ((a-b)*cos(i*π/(M+1))+(a+b))/2;
    end
    
    A = collocA(x);
    F = popF(x,α,β,f);

    C = A\F;

    U = zeros(M+2);
    U[1] = α;
    U[M+2] = β;
    for i ∈ 2:M+1
        U[i] = C ⋅ rowPop(x[i],M);
    end

    return U
end

# Main
trU(x) = exp(x) + (2-ℯ)*x - 1;

function comp(a,b,α,β,M)
    f(x) = exp(x);
    
    x = zeros(M+2);
    for i ∈ 0:M+1
        x[i+1] = ((a-b)*cos(i*π/(M+1))+(a+b))/2;
    end

    U = uSolve(a,b,α,β,f,M);

    a = plot(x,U-trU.(x));
    # a = plot!(x,trU.(x));
    display(a)
end

comp(0,1,0,1,20)
