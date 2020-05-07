### Package to handle all of the error calculations and plots of RBFT-FD
module Analysis

## Including used packages
using Plots
using LinearAlgebra: norm
using Arpack: eigs
gr()

# Public functions
export
    vecError,
    scalError,

    vectorPlot,
    nodePlot,
    errPlot,
    interPlot
# ----------------

## Error functions
# Relative vector error
function vecError(trues, calcs)
    N = size(trues,2);
    normErrs = zeros(N);
    mag = 0;
    for j ∈ 1:N
        magDiff = norm(trues[:,j]-calcs[:,j]);
        mag = max(norm(trues[:,j]),mag);
        normErrs[j] = magDiff;
    end
    normErrs /= mag;
    
    err = maximum(normErrs);
    
    return (err,normErrs)
end

# Relative scalar error
function scalError(trues, calcs)
    N = size(trues, 1);
    normErrs = zeros(N);
    mag = maximum(trues);
    for j ∈ 1:N
        magDiff = norm(trues[j]-calcs[j]);
        
        normErrs[j] = magDiff/mag;
    end
    err = maximum(normErrs);
    
    return (err,normErrs)
end

## Plot functions
# Plot of nodes and vectors
function vectorPlot(points, direction)
    N = size(points, 2);
    vecs = points + 0.1*direction;
    a = scatter(points[1,:],points[2,:],
                aspectratio = :equal,
                legend = false,
                show = false,
                markersize = 2,
                markerstrokealpha = 0,
                markeralpha = .9);
    for i ∈ 1:N
        plot!([points[1,i], vecs[1,i]],
              [points[2,i], vecs[2,i]],
              color = :red,
              legend = false,
              show = false,
              lw = 1)
    end
    return a
end

# Nodes plot
function nodePlot(set1)
    a = scatter(set1[1,:],set1[2,:],
                color = :blue,
                show = false,
                legend = false,
                aspectratio = :equal,
                markersize = 1,
                markerstrokealpha = 0)
    return a
end

# Error plot
function errPlot(nodes, errs)
    N = size(errs,1);

    # Scale errors logarithmically
    errs = log10.(errs);

    a = scatter(nodes[1,:], nodes[2,:],
                marker_z = errs,
                c = :viridis,
                cbarlims = :auto,
                colorbar = :right,
                aspectratio = :equal,
                legend = false,
                markersize = 2,
                markerstrokealpha = 0,
                markeralpha = .75,
                title = "Normal Error")
    return a
end

# RBF interpolation plot
# PHS
ϕ(x1,x2,m) = abs(x1-x2)^m;

# Interpolation function
function S(x1, x, λ, m)
    n = size(x,1);
    nn = size(λ,1) - n;
    nt = size(x1,1);
    
    s = zeros(nt);

    # Add RBF contribution
    for i ∈ 1:n, j ∈ 1:nt
        s[j] += λ[i]*ϕ(x1[j], x[i], m);
    end

    # Add polynomial contribution
    for i ∈ 0:nn-1, j ∈ 1:nt
        s[j] += λ[n+i+1]*x1[j]^i;
    end

    return s
end

function interPlot(nodes, λ, m)
    x = nodes[1,:];
    y = nodes[2,:];
    t = range(minimum(x), maximum(x), length = 100);
    s = S(t, x, λ, m);
    
    a = scatter(x,y,
                color = :blue,
                aspectratio = :equal,
                show = false,
                legend = false,
                markersize = 1.5,
                markerstrokealpha = 0)
    plot!(t,s,lw=1)

    return a
end

# Spectrum plot
function spectrum(A)
    # Find eigenvalues of A
    λ = eigs(A, maxiter = 1000, nev = 3000)[1];

    a = scatter(real(λ), imag(λ));

    return a
end

end # Module
