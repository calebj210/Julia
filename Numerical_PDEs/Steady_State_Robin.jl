### 1D Steady State Solver with Robin BC
using LinearAlgebra
using Plots

function popA(a,b,α,β,M)
    h = 1/(M+1);

    # Construct sub and super diagonal of A
    du = fill(1.0,M);
    dl = fill(1.0,M);
    du[1] = h;

    # Construct diagonal of A
    d = fill(-2.0,M+1);
    d[1] = -h;

    # Construct A
    A = Tridiagonal(dl, d, du);
    A /= h^2;
    
    return A
end

function popF(a,b,α,β,f,M)
    h = 1/(M+1);
    x = a:h:b-h;

    # Populate F
    F = f.(x);
    F[1] = α + F[1]*h/2;
    F[M+1] -= β/(h^2);
    
    return F
end

function solveODE()
    a = 0;
    α = 1;
    b = 1;
    β = -1;
    pts = 20*(2 .^ (1:7));
    f(x) = x^2;

    trU(x) = (x^4 + 12x - 25)/12;

    errs = zeros(size(pts,1));
    i = 1;
    for M ∈ pts
        A = popA(a,b,α,β,M);
        F = popF(a,b,α,β,f,M);

        U = A\F;
        
        h = 1/(M+1);
        x = a:h:b-h;

        err = U - trU.(x);
        errs[i] = norm(err)/(√(M));
        i += 1;
    end

    for i ∈ 2:size(pts,1)
        display(errs[i-1]/errs[i])
    end
end

solveODE()
