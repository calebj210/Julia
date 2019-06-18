### Finding Normals Using RBFs and Tensors

## Including used packages
using LinearAlgebra
using NearestNeighbors
using Plots
using BenchmarkTools
pyplot()

## Function definitions
# Unit Sphere
sph(u,v) = [cos(u)*sin(v), sin(u)*sin(v), cos(v)];

## Generate 3D Cartesian node set given u and v coordinate sets
function dist(u,v)
    N = size(u,1);
    M = size(v,1);

    nodes = ones(3,N*M);
    for i ∈ 1:N, j ∈ 1:M
        nodes[:, (j-1)*N+i] = sph(u[i], v[j]);
    end

    return nodes
end

## Generate uniform random node distribution on the unit sphere
function randSph(coords)
    x1 = coords[1];
    x2 = coords[2];

    tmp = x1^2 + x2^2;
    if tmp < 1
        x = 2*x1*sqrt(1-tmp);
        y = 2*x2*sqrt(1-tmp);
        z = 1- 2*tmp;
        return [x,y,z]
    else
        return randSph(rand(2)*2 .- 1)
    end
end


function randDist(N)
    nodes = zeros(3,N);

    # Create random distribution of coordinate pairs
    tmpNodes = rand(2,N)*2 .- 1;

    # Compute uniformly random unit sphere points using random
    # coordinate pairs
    for i ∈ 1:N
        nodes[:,i] = randSph(tmpNodes[:,i]);
    end

    return nodes
end

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
