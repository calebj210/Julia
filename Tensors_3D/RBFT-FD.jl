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
function appNorms(nodes, idx)
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

    return cents
end

# Positive Nth-D rotation function based on Givens rotations
function rotUp(nodes, normal)
    # Dimension of space
    D= size(nodes,1);
    cluster = copy(nodes);
    nml = copy(normal);
    
    # Begin rotations (under construction)
    sc = zeros(2,D-1)
    for i ∈ 1:D-1
        x = nml[i];
        y = nml[D];
        mag = norm([x,y])

        # Compute sine and cosine values
        s = x/mag;
        c = y/mag;
        
        # Store sine and cosine values for later
        sc[:,i] = [s,c];
        
        # Construct rotation matrix
        rot = [c -s;
               s c];
        
        nml[[i,D]] = [0,mag]
        cluster[[i,D],:] = rot*cluster[[i,D],:];
    end
    
    return (cluster,sc)
end

# Reverse rotation function given sine and cosines of previous rotations
function rotBack(vector, sc)
    D = size(vector,1);
    vec = copy(vector);
    
    for i ∈ D-1:-1:1
        s = sc[1,i];
        c = sc[2,i];
        
        # Construct rotation matrix
        rot = [c s;
               -s c];
        
        vec[[i,end]] = rot*vec[[i,end]];
    end

    return vec
end

## Polynomial functions
# Degree array generator
function polMat(vars, deg)
    if deg < 0
        return []
    else
        tmp::Array{Int64,1} = zeros(vars);
        array::Array{Array{Int64,1},1} = [copy(tmp)];
        idx = 1;
        i = 1;
        while true
            if tmp[1] == deg
                break
            end
            if sum(tmp) == deg
                tmp[idx] = 0;
                tmp[idx-1] += 1;
                idx -= 1;
                push!(array,copy(tmp));
            elseif idx != vars
                idx += 1;
            else
                tmp[idx] += 1;
                push!(array,copy(tmp));
            end
        end
        return reduce(hcat, array)
    end
end
