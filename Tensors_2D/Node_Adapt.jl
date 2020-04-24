### Preliminary adaptive node setup for curves
using ElasticArrays
using LinearAlgebra
include("/home/cajacobs/Documents/Julia/Tensors_2D/Packages/EvolvingCurves.jl/src/Evolving_Backend/Parameterization.jl")
include("/home/cajacobs/Documents/Julia/Tensors_2D/Packages/EvolvingCurves.jl/src/Evolving_Backend/RBF_Functions.jl")

## Function for tensor normals
function nodeAdapt!(nodes,n,cr=0.25,ca=0.75)
    # Find number of nodes
    N = size(nodes, 2);
    
    # Compute approximate normals
    appNorms = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFullOrd(nodes, n);
    
    normals = zeros(2,N);
    prox = zeros(n-1);
    λ = zeros(n);
    for i ∈ 1:N
        # Preallocate space for new nodes
        newNodes = [];
        
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Compute sorting permutation of  nodes from left to right
        p = sortperm(rot[1,:])

        totalDiff = 0;
        for j ∈ 1:n-1
            totalDiff += norm(rot[:,p[j+1]]-rot[:,p[j]]);
        end

        for j ∈ 1:n-1
            prox[j] = (norm(rot[:,p[j+1]]-rot[:,p[j]]))/totalDiff;
        end
        
        
    end

    nodes[1,1]=1000
    
    return normals
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
    display(nodes[1,1])

    nodeAdapt!(nodes, n)

    # display("Success")
    display(nodes[1,1])
end

main(100,10,3,2)
