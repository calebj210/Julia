### Preliminary adaptive node setup for curves
using ElasticArrays
using LinearAlgebra
using Analysis
include("/home/cajacobs/Documents/Julia/Tensors_2D/Packages/EvolvingCurves.jl/src/Evolving_Backend/Parameterization.jl")
include("/home/cajacobs/Documents/Julia/Tensors_2D/Packages/EvolvingCurves.jl/src/Evolving_Backend/RBF_Functions.jl")

## Function for adding and removing nodes
function nodeAdapt!(Nodes,n, m = 5, o = 3;cr=0.25,ca=0.75)
    # Reshape nodes into functional form
    nodes = reshape(copy(Nodes),2,:);

    # Find number of nodes
    N = size(nodes, 2);

    # Find nearest neighbors
    idx = knnFull(nodes, n);
    
    # Compute approximate normals
    appNorms = approxNormals(nodes, idx);
    appNorms = orVecs(nodes, appNorms, idx);
    if appNorms[:,1]⋅[1,0] < 0
        appNorms *= -1;
    end

    # C = vectorPlot(nodes, appNorms)
    # display(C)
    # sleep(0.5)

    prox = zeros(n-1);
    newNodes = [];
    for i ∈ 1:N
        # Find the angle between the normal vector and the ĵ direction
        θ = findAngle(appNorms[:,i]);

        # Center nodes in cluster around the node of interest
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes about the center node by θ rads
        rot = rotate(cent, θ);

        # Compute sorting permutation of  nodes from left to right
        p = sortperm(rot[1,:])

        # Compute proxy at each subinterval
        totalDiff = 0;
        for j ∈ 1:n-1
            totalDiff += norm(rot[:,p[j+1]]-rot[:,p[j]]);
        end

        for j ∈ 1:n-1
            prox[j] = (norm(rot[:,p[j+1]]-rot[:,p[j]]))/totalDiff;
        end

        # Find interpolation weights
        λ = interpolate(rot, m, o);

        # Find the x1 coordinate of the node halfway between 0 and the
        # next node to the right
        x1 = 0;
        for j = 1:n-1
            if rot[1,p[j]] == 0
                x1 = rot[1,p[j+1]]/2;
                break
            end
        end

        # Compute x2 (the height) at x1 given our interpolation
        x2 = S(x1, rot[1,:], λ, m);

        # Compute added node
        x = [x1; x2];
        x = rotate(x, -θ);
        x += nodes[:,i];
        
        # Append new node to our full node set
        append!(Nodes, x);
    end
end

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
    Nodes = nodes[:];

    idx = knnFull(nodes, n)

    # Compute approximate normals
    appNorms = approxNormals(nodes, idx);
    
    # for i = 1:N
    #     # Find the angle between the normal vector and the ĵ direction
    #     θ = findAngle(appNorms[:,i]);

    #     # Center nodes in cluster around the node of interest
    #     cent = center(nodes[:,idx[i]]);

    #     # Rotate nodes about the center node by θ rads
    #     rot = rotate(cent, θ);

    #     # Find interpolation weights
    #     λ = interpolate(rot, m, o);
        
    #     C = interPlot(rot, λ, m)
    #     display(C)
    #     sleep(0.1)
    # end

    A = nodePlot(nodes);
    
    nodeAdapt!(Nodes, n, m, o);
    
    nodes = reshape(Nodes, (2,:));
    
    errs = []
    for i∈1:size(nodes,2)
        push!(errs,abs(norm(nodes[:,i])-1))
    end

    # B = nodePlot(reshape(Nodes, 2, :));

    # display(A)
    # sleep(1)
    # display(B)

    display(errPlot(nodes, errs))
end

main(11,5,5,2)
