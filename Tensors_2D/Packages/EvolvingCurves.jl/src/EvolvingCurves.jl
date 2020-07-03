module EvolvingCurves
### RBFT-FD for Evolving Curves
## Functions to export
export
    findNormals,
    computeκ,
    ∇∇,
    constructLBO,
    constructBHO
# -----------------------------

## Including used packages
using LinearAlgebra
include("Evolving_Backend/Parameterization.jl")
include("Evolving_Backend/RBF_Functions.jl")

## Function for tensor normals
function findNormals(nodes, n, m, o)
    # Find number of nodes
    N = size(nodes, 2);
   
    # Find nearest neighbors
    idx = knnFull(nodes, n);

    # Compute approximate normals
    appNorms = approxNormals(nodes, idx);

    # Begin computing normals at each node
    normals = zeros(2,N);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Compute RBF weights
        λ = interpolate(rot, m, o);

        # Compute derivative at S1 = 0
        f_s = S_x(rot[1,:], λ, m);
        
        # Compute length element
        s = 1 + f_s^2;
        L = sqrt(s);
        
        # Compute normal at S1 = 0
        tmp = [-f_s/L, 1/L];
        
        # Rotate local normal to proper angle
        normals[:,i] = rotate(tmp, -θ);
    end

    normals = orVecs(nodes, normals, idx)
    
    return normals
end

## Curvature Computation
function computeκ(nodes, n, m, o)
    # Find number of nodes
    N = size(nodes, 2);
    
    # Compute approximate normals
    appNorms = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFull(nodes, n);

    # Compute curvature at each node
    κ = zeros(N,1)
    for i ∈ 1:N
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Create collocation matrix
        A = collocM(rot[1,:], m, o);
        
        # Create local height vector with helping terms
        f = [rot[2,:]; zeros(o+1)];
        
        # Compute RBF weights of local surface
        λs = A\f;

        # Compute surface derivatives at each local node S = 0
        S_s = S_x(rot[1,:], λs, m);
        S_ss = S_xx(rot[1,:], λs, m);
        
        # Initialize linear operator vector
        κ[i] = abs(S_ss)/((1+S_s^2)^(3/2));
    end
    
    return κ
end

## Laplace-Beltrami Operator
# This is for directly computing the LBO
function ∇∇(nodes, F, n, m, o)
    # Find number of nodes
    N = size(nodes, 2);
    
    # Compute approximate normals
    appNorms = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFull(nodes, n);

    # Allocate space for Laplace-Beltrami values
    vals = zeros(N);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Compute RBF weights of local surface
        λs = interpolate(rot, m, o);

        # Compute RBF weights of local function
        λf = interpolate([rot[1,:]'; F[idx[i]]'], m, o);

        # Compute derivatives at S1 = 0
        S_s = S_x(rot[1,:], λs, m);
        S_ss = S_xx(rot[1,:], λs, m);
        S_f = S_x(rot[1,:], λf, m);
        S_ff = S_xx(rot[1,:], λf, m);
        
        # Compute length element
        s = 1 + S_s^2;
        
        # Compute ∇∇F
        vals[i] = s^(-1)*S_ff - s^(-2)*S_s*S_ss*S_f;
    end
    
    return vals
end

## Discretizations
# Populate Laplace-Beltrami discretization matrix
function constructLBO(nodes, n, m, o)
    # Find number of nodes
    N = size(nodes, 2);
    
    # Compute approximate normals
    appNorms = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFull(nodes, n);
    
    cent = zeros(2,n);
    rot = cent;
    D = spzeros(N,N);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Create collocation matrix
        A = collocM(rot[1,:], m, o);
        
        # Create local height vector with helping terms
        f = [rot[2,:]; zeros(o+1)];
        
        # Compute RBF weights of local surface
        λs = A\f;

        # Compute derivatives at each local node S = 0
        S_s = S_x(rot[1,:], λs, m);
        S_ss = S_xx(rot[1,:], λs, m);
        
        # Compute length element
        s = 1 + S_s^2;
        
        # Initialize linear operator vector
        Lϕ = zeros(n+o+1);

        # Populate non-polynomial vector components
        for i ∈ 1:n
            xi = rot[1,i];

            # Laplace-Beltrami operator
            Lϕ[i] = s^(-1)*ϕ_xx(0, xi, m) - s^(-2)*S_s*S_ss*ϕ_x(0, xi, m);
        end
        
        # Populate polynomial vector components
        if o > 0
            Lϕ[n+2] = -s^(-2)*S_s*S_ss;
            if o > 1
                Lϕ[n+3] = 2*s^(-1);
            end
        end

        # Compute local weights
        w = A\Lϕ;
        
        # Populate D with local weights
        D[i,idx[i]] = w[1:n];
    end
    
    return D
end

# Populate Biharmonic Operator discretization matrix
function constructBHO(nodes, n, m, o)
    # Find number of nodes
    N = size(nodes, 2);
    
    # Compute approximate normals
    appNorms = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFull(nodes, n);
    
    cent = zeros(2,n);
    rot = cent;
    D = spzeros(N,N);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Create collocation matrix
        A = collocM(rot[1,:], m, o);
        
        # Create local height vector with helping terms
        f = [rot[2,:]; zeros(o+1)];
        
        # Compute RBF weights of local surface
        λs = A\f;

        # Compute derivatives at each local node S = 0
        S_s = S_x(rot[1,:], λs, m);
        S_ss = S_xx(rot[1,:], λs, m);
        if typ == 2
            S_sss = S_xxx(rot[1,:], λs, m);
            S_ssss = S_xxxx(rot[1,:], λs, m);
        end
        
        # Compute length element
        s = 1 + S_s^2;
        
        # Initialize linear operator vector
        Lϕ = zeros(n+o+1);

        # Populate non-polynomial part of linear operator vector
        for i ∈ 1:n
            xi = rot[1,i];

            
            # Biharmonic operator
            Lϕ[i] =
                (13*s^(-4)*S_s*S_ss^3 - 28*s^(-5)*S_s^3*S_ss^3 -
                 3*s^(-3)*S_ss*S_sss + 13*s^(-4)*S_s^2*S_ss*S_sss -
                 s^(-3)*S_s*S_ssss)*
            ϕ_x(0, xi, m) + 
                (19*s^(-4)*S_s^2*S_ss^2 - 4*s^(-3)*S_ss^2 -
                 4*s^(-3)*S_s*S_sss)*
            ϕ_xx(0, xi, m) +
                (-6*s^(-3)*S_s*S_ss)*
            ϕ_xxx(0, xi, m) +
                (s^(-2))*
            ϕ_xxxx(0, xi, m);
        end

        # Populate polynomial part of linear operator vector
        if o > 0
            Lϕ[n+2] =
                (13*s^(-4)*S_s*S_ss^3 - 28*s^(-5)*S_s^3*S_ss^3 -
                 3*s^(-3)*S_ss*S_sss + 13*s^(-4)*S_s^2*S_ss*S_sss -
                 s^(-3)*S_s*S_ssss);
            if o > 1
                Lϕ[n+3] =
                    2*(19*s^(-4)*S_s^2*S_ss^2 -
                       4*s^(-3)*S_ss^2 -
                       4*s^(-3)*S_s*S_sss);
                if o > 2
                    Lϕ[n+4] =
                        6*(-6*s^(-3)*S_s*S_ss);
                    if o > 3
                        Lϕ[n+5] =
                            24*s^(-2);
                    end
                end
            end
        end
        
        # Compute local weights
        w = A\Lϕ;
        
        # Populate D with local weights
        D[i,idx[i]] = w[1:n];
    end
    
    return D
end

end # module
