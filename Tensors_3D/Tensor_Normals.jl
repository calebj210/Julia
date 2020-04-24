### Finding Normals Using RBFs and Tensors

## Including used packages and files
using LinearAlgebra
using NearestNeighbors
using Plots
using BenchmarkTools
include("Nodes.jl")
pyplot()

## Full function
# For exact control of each coordinate
function comp(N, M)
    u = range(0, 2*π-2*π/N, length = N);
    v = range(0+π/(2*M), π-π/(2*M), length = M);

    nodes =  dist(u, v);

    plot(nodes[1,:], nodes[2,:], nodes[3,:], l=:scatter3d)
end

# Uses roughly NN number of nodes for our surface
function comp(NN)
    N = round(Int,sqrt(NN));
    u = range(0, 2*π-2*π/N, length = N);
    v = range(0+π/(2*N), π-π/(2*N), length = N);
    
    nodes =  dist(u, v);
    
    plot(nodes[1,:], nodes[2,:], nodes[3,:], l=:scatter3d)
end

# Uses N random nodes on the unit sphere
function rndComp(N)
    nodes = randDist(N);

    scatter(nodes[1,:], nodes[2,:], nodes[3,:])
end

rndComp(500)
