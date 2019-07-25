### RBFT-FD functions

## Including used packages
using NearestNeighbors
using LinearAlgebra

### Approximate normal functions definitions
# Inverse iteration routine to find eigenvector assoiciated with smallest λ
function invIt(A; ϵ = 10^(-5), maxIts = 10)
    # If A is singular bypass inverse iteration and return the nullspace instead
    if det(A) != 0
        B = lu(A);
        x = rand(size(A,1))
        y1 = x;
        for i ∈ 1:maxIts
            y0 = y1;
            y1 = B\x;
            x = y1/norm(y1, Inf);
            
            if norm(y0-y1, Inf) < ϵ || norm(y0+y1, Inf) < ϵ
                return x
            end
        end
    else
        return nullspace(A)[:,1]
    end

    return x
end

# Generate covariance given a cluster xⱼ
function covar(nodes)
    # Number of nodes
    N = size(nodes,2);
    
    # Find average vector
    x̄ = nodes[:,1]
    for i ∈ 2:N
        x̄ += nodes[:,2];
    end
    x̄ = x̄ / N;

    # Compute covariance
    A = zeros(size(nodes,1),size(nodes,1))
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
function appNorms(nodes, idx; n::Int = 3)
    # Number of nodes
    N = size(nodes,2);

    # Compute approximate normal by returning the eigenvector
    # associated with the smallest λ of the covariance matrix
    nmls = zeros(size(nodes,1), N);
    for i ∈ 1:N
        c = covar(nodes[:,idx[i]])
        nmls[:,i] = invIt(c);
    end

    # Orientation (under construction)

    return nmls
end


## Node orienting functions
# Node centering function given where offset is the first node in the array
function center(nodes)
    cents = nodes .- nodes[:,1];

    return nodes
end

# Positive Nth-D rotation function based on Givens rotations
function rotUp(nodes, nml, n::Int = 3)
    # Dimension of space
    D= size(nodes,1);
    
    # Begin rotations (under construction)
    sc = zeros(2,N-1)
    for i ∈ 1:D-1
        sc[:,i] = [nml]
    end

    return nodes
end

# Reverse rotation function given sine and cosines of previous rotations
function rotBack(nodes, sc)
    
end
