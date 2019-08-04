### Calculation functions
## Including used packages and files
using Plots
include("RBFT-FD.jl")
include("Nodes.jl")
pyplot()

## Function for finding 3D surface normals
function find3DNormals(nodes, n, m, deg; idx = [])
    # Number of nodes
    N = size(nodes, 2);

    # Determine if KNN is needed for finding node clusters
    if idx == []
        idx = knnFull(nodes, n);
    end

    # Find approximate normals
    appNrms = appNorms(nodes, idx);

    # Create polynomial power matrix;
    poly = polyMat(2, deg);
    
    # Compute better normals
    normals = zeros(3,N);
    for i ∈ 1:N
        # Grab local cluster
        nds = copy(nodes[:,idx[i]]);

        # Center cluster around cartesian origin
        center!(nds);

        # Rotate nodes to have approximate normal point in the positve z direction
        sc = rotUp!(nds, appNrms[:,i]);

        # Find local interpolant weights
        λ = findλ(nds[1:2, :], nds[3,:], poly, m = m);

        # Compute local derivatives
        S_x = S_xi(nds[1:2, :], [0,0], λ, 1, poly, m = m);
        S_y = S_xi(nds[1:2, :], [0,0], λ, 2, poly, m = m);
        
        # Compute better normal
        tmpN = [-S_x, -S_y, 1];
        tmpN = tmpN/norm(tmpN);
        normals[:,i] = rotBack(tmpN,sc);
    end

    return normals
end

## Test for find3DNormals
function comp()
    nodes = randDist(1000);

    nmls = find3DNormals(nodes, 5, 3, 1);

    a = scatter(nodes[1,:], nodes[2,:], nodes[3,:]);
    for i ∈ 1:size(nodes,2)
        a = plot!([nodes[1,i],nodes[1,i]+.5*nmls[1,i]],
                  [nodes[2,i],nodes[2,i]+.5*nmls[2,i]],
                  [nodes[3,i],nodes[3,i]+.5*nmls[3,i]],
                  linecolor = :red,
                  legend = false);
    end
    
    display(a)
end

comp()
