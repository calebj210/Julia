### Parameterization Essentials

## Including used packages
using NearestNeighbors

## Functions for generating approximate normals using a sorted node
## distribution
# function for finding approximate normals
function approxNormals(nodes)
    # Number of nodes
    N = size(nodes,2);
    nmls = zeros(2,N);
    
    # First Node
    tmp = nodes[:,2] - nodes[:,N];
    tmp = rot90(tmp);
    nmls[:,1] = tmp;

    # Middle Nodes
    for i ∈ 2:N-1
        tmp = nodes[:,i+1] - nodes[:, i-1];
        tmp = rot90(tmp);
        nmls[:,i] = tmp;
    end

    # Last Node
    tmp = nodes[:,1] - nodes[:,N-1];
    tmp = rot90(tmp);
    nmls[:,N] = tmp;

    return nmls
end

# Function that rotates vectors -π/2
function rot90(vec)
    tmp = zeros(2);
    tmp[1] = vec[2];
    tmp[2] = -vec[1];
    
    return tmp
end

## Find indices for k nearest neighbors
# KNN for an unordered data set
function knnFull(nodes, n)
    kdtree = KDTree(nodes);
    idx, tmp = knn(kdtree, nodes, n, true);
    
    return idx
end

# KNN for ordered data set
function knnFullOrd(nodes, n)::Array{Array{Int,1},1}
    N = size(nodes,2);
    idx = [];
    for i ∈ 1:N
        push!(idx,zeros(n))
    end
    
    # Put nth node as nth first index
    for i ∈ 1:N
        idx[i][1] = i;
    end

    # Populate rest of indices
    for i ∈ 1:N, j ∈ 2:n
        if j%2 == 0
            k = i + floor(j/2);
            idx[i][j] = k>N ? k - N : k;
        else
            k = i - floor(j/2);
            idx[i][j] = k<=0 ? N+k : k;
        end 
    end
    
    return idx
end

## enable ordered knn search
# knnFull(nodes, n) = knnFullOrd(nodes, n);

## Functions for setting up RBF envrionment
# Center main node
function center(nodes)
    nodes = nodes .- nodes[:,1];

    return nodes
end

function unCenter(nodes,center)
    node 
    
    return nodes
end

# Finds angle between center vector and ĵ
sgn(x) = x >= 0 ? 1 : -1;
function findAngle(node)
    b = node[2]/norm(node);
    θ = sgn(node[1])acos(b);
    
    return θ
end

# Rotate nodes about origin
function rotate(nodes, θ)
    # Compute rotation matrix
    rotor = [cos(θ) -sin(θ);
             sin(θ) cos(θ)];
    
    # Rotate all vectors
    nodes = rotor*nodes;
    
    return nodes
end
