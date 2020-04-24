### Computing Normals

## Including used packages
using EvolvingCurves: findNormals
using Analysis

## Unit circle
x(t) = cos.(t);
y(t) = sin.(t);

## True normals
truX(t) = cos.(t);
truY(t) = sin.(t);

## Start
function main(N=20, n=5, m=3, o=0)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N);

    # Creating our node distribution
    nodes = [x(t)';y(t)'];

    # Compute normals
    norms = findNormals(nodes, n, m, o);

    # Computing true normals
    truNorms = [truX(t)';truY(t)'];

    # Plot of results
    plotA = vectorPlot(nodes, norms);

    # Error calculation
    error,tmp = vecError(truNorms, norms);

    display(plotA)

    return error
end

main(50,5,3,2)
