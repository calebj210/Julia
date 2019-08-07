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
    # N = 23;
    # M = 12;
    # θ = range(0, 2*π*(1-1/N), length = N);
    # ϕ = range(π/4, 3*π/4, length = M);
    # nodes = dist(θ,ϕ);
    nodes = randDist(15000);

    deg = 2;
    
    tmp = size(polyMat(2,deg),2);
    nmls = find3DNormals(nodes,  2*tmp, 3, deg);
    appnorms = appNorms(nodes, knnFull(nodes,2*tmp));
    for i ∈ 1:size(nodes,2)
        appnorms[:,i] /= norm(appnorms[:,i]);
    end
    for i ∈ 1:size(nodes,2)
        if nmls[:,i] ⋅ nodes[:,i] < 0
            nmls[:,i] *= -1;
        end
        if appnorms[:,i] ⋅ nodes[:,i] < 0
            appnorms[:,i] *= -1;
        end
    end

    appnorms -= nodes;
    nmls -= nodes;
    nms1 = zeros(size(nodes,2));
    nms2 = zeros(size(nodes,2));
    for i ∈ 1:size(nodes,2);
        nms1[i] = norm(nmls[:,i]);
        nms2[i] = norm(appnorms[:,i]);
    end

    display(maximum(nms1))
    display(maximum(nms2))
    
    # a = scatter(nodes[1,:], nodes[2,:], nodes[3,:]);
    # for i ∈ 1:size(nodes,2)
    #     a = plot!([nodes[1,i],nodes[1,i]+nmls[1,i]],
    #               [nodes[2,i],nodes[2,i]+nmls[2,i]],
    #               [nodes[3,i],nodes[3,i]+nmls[3,i]],
    #               linecolor = :green,
    #               legend = false);
    # end
    # for i ∈ 1:size(nodes,2)
    #     a = plot!([nodes[1,i],nodes[1,i]+appnorms[1,i]],
    #               [nodes[2,i],nodes[2,i]+appnorms[2,i]],
    #               [nodes[3,i],nodes[3,i]+appnorms[3,i]],
    #               linecolor = :red,
    #               legend = false);
    # end
    
    # display(a)
end

comp()
