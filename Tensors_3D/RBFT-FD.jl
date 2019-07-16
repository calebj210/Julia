### RBFT-FD functions

## Including used packages
using NearestNeighbors

### Approximate normal functions definitions
# Generate covariance given a cluster xⱼ
function covar(nodes)
    N = size(nodes,2);
    
    # Find average vector
    x̄ = nodes[:,1]
    for i ∈ 2:N
        x̄ += nodes[:,2];
    end
    x̄ = x̄ / N;

    # Compute covariance
    A = zeros(N,N)
    for i ∈ 1:N
        diff = nodes[:,i] - x̄;

        A += diff*diff';
    end
    
    return A
end

# KNN search for unordered data set
function knnFull(nodes, n)
    kdtree = KDTree(nodes);
    idx, tmp = knn(kdtree, nodes, n, true);
    
    return idx
end

# Approximate normals using associated covariance
function appNorms(nodes; n::Int = 3)
    N = size(nodes,2);

    # Find nearest neighbors
    idx = knnFull(nodes, n);
    
    nmls = zeros(size(nodes,1), N);
    for i ∈ 1:N
        c = covar(nodes[:,idx[i]])
        nmls[:,i] = eigvecs(c, eigmin(c));
    end
end
