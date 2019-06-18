### Finding Normals Using RBFs and Tensors

## Including used packages
using LinearAlgebra
using SparseArrays
using NearestNeighbors
using Plots
using BenchmarkTools
gr()


## Function definitions
# Circle
circX(t) = cos.(t);
circY(t) = sin.(t);

# True circle normals
circNormsX(t) = cos.(t);
circNormsY(t) = sin.(t);

# Pi curve
piX(t) = ((17/31)*sin.((235/57)+(-32)*t)+(19/17)*sin.((192/55)+(-30)*t)+(47/32)*sin.((69/25)+(-29)*t)+(35/26)*sin.((75/34)+(-27)*t)+(6/31)*sin.((23/10)+(-26)*t)+(35/43)*sin.((10/33)+(-25)*t)+(126/43)*sin.((421/158)+(-24)*t)+(143/57)*sin.((35/22)+(-22)*t)+(106/27)*sin.((84/29)+(-21)*t)+(88/25)*sin.((23/27)+(-20)*t)+(74/27)*sin.((53/22)+(-19)*t)+(44/53)*sin.((117/25)+(-18)*t)+(126/25)*sin.((88/49)+(-17)*t)+(79/11)*sin.((43/26)+(-16)*t)+(43/12)*sin.((41/17)+(-15)*t)+(47/27)*sin.((244/81)+(-14)*t)+(8/5)*sin.((79/19)+(-13)*t)+(373/46)*sin.((109/38)+(-12)*t)+(1200/31)*sin.((133/74)+(-11)*t)+(67/24)*sin.((157/61)+(-10)*t)+(583/28)*sin.((13/8)+(-8)*t)+(772/35)*sin.((59/16)+(-7)*t)+(3705/46)*sin.((117/50)+(-6)*t)+(862/13)*sin.((19/8)+(-5)*t)+(6555/34)*sin.((157/78)+(-3)*t)+(6949/13)*sin.((83/27)+(-1)*t)+(-6805/54)*sin.((1/145)+2*t)+(-5207/37)*sin.((49/74)+4*t)+(-1811/58)*sin.((55/43)+9*t)+(-63/20)*sin.((2/23)+23*t)+(-266/177)*sin.((13/18)+28*t)+(-2/21)*sin.((7/16)+31*t))/1000;

piY(t) =
    ((70/37)*sin.((65/32)+(-32)*t)+(11/12)*sin.((98/41)+(-31)*t)+(26/29)*sin.((35/12)+(-30)*t)+(54/41)*sin.((18/7)+(-29)*t)+(177/71)*sin.((51/19)+(-27)*t)+(52/33)*sin.((133/52)+(-26)*t)+(59/34)*sin.((125/33)+(-26)*t)+(98/29)*sin.((18/11)+(-25)*t)+(302/75)*sin.((59/22)+(-24)*t)+(104/9)*sin.((118/45)+(-22)*t)+(52/33)*sin.((133/52)+(-21)*t)+(37/45)*sin.((61/14)+(-20)*t)+(143/46)*sin.((144/41)+(-19)*t)+(254/47)*sin.((19/52)+(-18)*t)+(246/35)*sin.((92/25)+(-17)*t)+(722/111)*sin.((176/67)+(-16)*t)+(136/23)*sin.((3/19)+(-15)*t)+(273/25)*sin.((32/21)+(-13)*t)+(229/33)*sin.((117/28)+(-12)*t)+(19/4)*sin.((43/11)+(-11)*t)+(135/8)*sin.((23/10)+(-10)*t)+(205/6)*sin.((33/23)+(-8)*t)+(679/45)*sin.((55/12)+(-7)*t)+(101/8)*sin.((11/12)+(-6)*t)+(2760/59)*sin.((40/11)+(-5)*t)+(1207/18)*sin.((21/23)+(-4)*t)+(8566/27)*sin.((39/28)+(-3)*t)+(12334/29)*sin.((47/37)+(-2)*t)+(15410/39)*sin.((185/41)+(-1)*t)+(-596/17)*sin.((3/26)+9*t)+(-247/28)*sin.((25/21)+14*t)+(-458/131)*sin.((21/37)+23*t)+(-41/36)*sin.((7/8)+28*t))/1000;

# True pi normals
piNormsX(t) = (1/2000)*((-2240/37)*cos((65/32)+(-32)*t)+(-341/12)*cos((98/41)+(-31)*t)+(-780/29)*cos((35/12)+(-30)*t)+(-1566/41)*cos((18/7)+(-29)*t)+(-4779/71)*cos((51/19)+(-27)*t)+(-1352/33)*cos((133/52)+(-26)*t)+(-767/17)*cos((125/33)+(-26)*t)+(-2450/29)*cos((18/11)+(-25)*t)+(-2416/25)*cos((59/22)+(-24)*t)+(-2288/9)*cos((118/45)+(-22)*t)+(-364/11)*cos((133/52)+(-21)*t)+(-148/9)*cos((61/14)+(-20)*t)+(-2717/46)*cos((144/41)+(-19)*t)+(-4572/47)*cos((19/52)+(-18)*t)+(-4182/35)*cos((92/25)+(-17)*t)+(-11552/111)*cos((176/67)+(-16)*t)+(-2040/23)*cos((3/19)+(-15)*t)+(-3549/25)*cos((32/21)+(-13)*t)+(-916/11)*cos((117/28)+(-12)*t)+(-209/4)*cos((43/11)+(-11)*t)+(-675/4)*cos((23/10)+(-10)*t)+(-820/3)*cos((33/23)+(-8)*t)+(-4753/45)*cos((55/12)+(-7)*t)+(-303/4)*cos((11/12)+(-6)*t)+(-13800/59)*cos((40/11)+(-5)*t)+(-2414/9)*cos((21/23)+(-4)*t)+(-8566/9)*cos((39/28)+(-3)*t)+(-24668/29)*cos((47/37)+(-2)*t)+(-15410/39)*cos((185/41)+(-1)*t)+(-5364/17)*cos((3/26)+9*t)+(-247/2)*cos((25/21)+14*t)+(-10534/131)*cos((21/37)+23*t)+(-287/9)*cos((7/8)+28*t))*((1/1000000)*((-2240/37)*cos((65/32)+(-32)*t)+(-341/12)*cos((98/41)+(-31)*t)+(-780/29)*cos((35/12)+(-30)*t)+(-1566/41)*cos((18/7)+(-29)*t)+(-4779/71)*cos((51/19)+(-27)*t)+(-1352/33)*cos((133/52)+(-26)*t)+(-767/17)*cos((125/33)+(-26)*t)+(-2450/29)*cos((18/11)+(-25)*t)+(-2416/25)*cos((59/22)+(-24)*t)+(-2288/9)*cos((118/45)+(-22)*t)+(-364/11)*cos((133/52)+(-21)*t)+(-148/9)*cos((61/14)+(-20)*t)+(-2717/46)*cos((144/41)+(-19)*t)+(-4572/47)*cos((19/52)+(-18)*t)+(-4182/35)*cos((92/25)+(-17)*t)+(-11552/111)*cos((176/67)+(-16)*t)+(-2040/23)*cos((3/19)+(-15)*t)+(-3549/25)*cos((32/21)+(-13)*t)+(-916/11)*cos((117/28)+(-12)*t)+(-209/4)*cos((43/11)+(-11)*t)+(-675/4)*cos((23/10)+(-10)*t)+(-820/3)*cos((33/23)+(-8)*t)+(-4753/45)*cos((55/12)+(-7)*t)+(-303/4)*cos((11/12)+(-6)*t)+(-13800/59)*cos((40/11)+(-5)*t)+(-2414/9)*cos((21/23)+(-4)*t)+(-8566/9)*cos((39/28)+(-3)*t)+(-24668/29)*cos((47/37)+(-2)*t)+(-15410/39)*cos((185/41)+(-1)*t)+(-5364/17)*cos((3/26)+9*t)+(-247/2)*cos((25/21)+14*t)+(-10534/131)*cos((21/37)+23*t)+(-287/9)*cos((7/8)+28*t))^2+(1/1000000)*((-544/31)*cos((235/57)+(-32)*t)+(-570/17)*cos((192/55)+(-30)*t)+(-1363/32)*cos((69/25)+(-29)*t)+(-945/26)*cos((75/34)+(-27)*t)+(-156/31)*cos((23/10)+(-26)*t)+(-875/43)*cos((10/33)+(-25)*t)+(-3024/43)*cos((421/158)+(-24)*t)+(-3146/57)*cos((35/22)+(-22)*t)+(-742/9)*cos((84/29)+(-21)*t)+(-352/5)*cos((23/27)+(-20)*t)+(-1406/27)*cos((53/22)+(-19)*t)+(-792/53)*cos((117/25)+(-18)*t)+(-2142/25)*cos((88/49)+(-17)*t)+(-1264/11)*cos((43/26)+(-16)*t)+(-215/4)*cos((41/17)+(-15)*t)+(-658/27)*cos((244/81)+(-14)*t)+(-104/5)*cos((79/19)+(-13)*t)+(-2238/23)*cos((109/38)+(-12)*t)+(-13200/31)*cos((133/74)+(-11)*t)+(-335/12)*cos((157/61)+(-10)*t)+(-1166/7)*cos((13/8)+(-8)*t)+(-772/5)*cos((59/16)+(-7)*t)+(-11115/23)*cos((117/50)+(-6)*t)+(-4310/13)*cos((19/8)+(-5)*t)+(-19665/34)*cos((157/78)+(-3)*t)+(-6949/13)*cos((83/27)+(-1)*t)+(-6805/27)*cos((1/145)+2*t)+(-20828/37)*cos((49/74)+4*t)+(-16299/58)*cos((55/43)+9*t)+(-1449/20)*cos((2/23)+23*t)+(-7448/177)*cos((13/18)+28*t)+(-62/21)*cos((7/16)+31*t))^2)^(-1);

piNormsY(t) = (-1/2000)*((1/1000000)*((-2240/37)*cos((65/32)+(-32)*t)+(-341/12)*cos((98/41)+(-31)*t)+(-780/29)*cos((35/12)+(-30)*t)+(-1566/41)*cos((18/7)+(-29)*t)+(-4779/71)*cos((51/19)+(-27)*t)+(-1352/33)*cos((133/52)+(-26)*t)+(-767/17)*cos((125/33)+(-26)*t)+(-2450/29)*cos((18/11)+(-25)*t)+(-2416/25)*cos((59/22)+(-24)*t)+(-2288/9)*cos((118/45)+(-22)*t)+(-364/11)*cos((133/52)+(-21)*t)+(-148/9)*cos((61/14)+(-20)*t)+(-2717/46)*cos((144/41)+(-19)*t)+(-4572/47)*cos((19/52)+(-18)*t)+(-4182/35)*cos((92/25)+(-17)*t)+(-11552/111)*cos((176/67)+(-16)*t)+(-2040/23)*cos((3/19)+(-15)*t)+(-3549/25)*cos((32/21)+(-13)*t)+(-916/11)*cos((117/28)+(-12)*t)+(-209/4)*cos((43/11)+(-11)*t)+(-675/4)*cos((23/10)+(-10)*t)+(-820/3)*cos((33/23)+(-8)*t)+(-4753/45)*cos((55/12)+(-7)*t)+(-303/4)*cos((11/12)+(-6)*t)+(-13800/59)*cos((40/11)+(-5)*t)+(-2414/9)*cos((21/23)+(-4)*t)+(-8566/9)*cos((39/28)+(-3)*t)+(-24668/29)*cos((47/37)+(-2)*t)+(-15410/39)*cos((185/41)+(-1)*t)+(-5364/17)*cos((3/26)+9*t)+(-247/2)*cos((25/21)+14*t)+(-10534/131)*cos((21/37)+23*t)+(-287/9)*cos((7/8)+28*t))^2+(1/1000000)*((-544/31)*cos((235/57)+(-32)*t)+(-570/17)*cos((192/55)+(-30)*t)+(-1363/32)*cos((69/25)+(-29)*t)+(-945/26)*cos((75/34)+(-27)*t)+(-156/31)*cos((23/10)+(-26)*t)+(-875/43)*cos((10/33)+(-25)*t)+(-3024/43)*cos((421/158)+(-24)*t)+(-3146/57)*cos((35/22)+(-22)*t)+(-742/9)*cos((84/29)+(-21)*t)+(-352/5)*cos((23/27)+(-20)*t)+(-1406/27)*cos((53/22)+(-19)*t)+(-792/53)*cos((117/25)+(-18)*t)+(-2142/25)*cos((88/49)+(-17)*t)+(-1264/11)*cos((43/26)+(-16)*t)+(-215/4)*cos((41/17)+(-15)*t)+(-658/27)*cos((244/81)+(-14)*t)+(-104/5)*cos((79/19)+(-13)*t)+(-2238/23)*cos((109/38)+(-12)*t)+(-13200/31)*cos((133/74)+(-11)*t)+(-335/12)*cos((157/61)+(-10)*t)+(-1166/7)*cos((13/8)+(-8)*t)+(-772/5)*cos((59/16)+(-7)*t)+(-11115/23)*cos((117/50)+(-6)*t)+(-4310/13)*cos((19/8)+(-5)*t)+(-19665/34)*cos((157/78)+(-3)*t)+(-6949/13)*cos((83/27)+(-1)*t)+(-6805/27)*cos((1/145)+2*t)+(-20828/37)*cos((49/74)+4*t)+(-16299/58)*cos((55/43)+9*t)+(-1449/20)*cos((2/23)+23*t)+(-7448/177)*cos((13/18)+28*t)+(-62/21)*cos((7/16)+31*t))^2)^(-1)*((-544/31)*cos((235/57)+(-32)*t)+(-570/17)*cos((192/55)+(-30)*t)+(-1363/32)*cos((69/25)+(-29)*t)+(-945/26)*cos((75/34)+(-27)*t)+(-156/31)*cos((23/10)+(-26)*t)+(-875/43)*cos((10/33)+(-25)*t)+(-3024/43)*cos((421/158)+(-24)*t)+(-3146/57)*cos((35/22)+(-22)*t)+(-742/9)*cos((84/29)+(-21)*t)+(-352/5)*cos((23/27)+(-20)*t)+(-1406/27)*cos((53/22)+(-19)*t)+(-792/53)*cos((117/25)+(-18)*t)+(-2142/25)*cos((88/49)+(-17)*t)+(-1264/11)*cos((43/26)+(-16)*t)+(-215/4)*cos((41/17)+(-15)*t)+(-658/27)*cos((244/81)+(-14)*t)+(-104/5)*cos((79/19)+(-13)*t)+(-2238/23)*cos((109/38)+(-12)*t)+(-13200/31)*cos((133/74)+(-11)*t)+(-335/12)*cos((157/61)+(-10)*t)+(-1166/7)*cos((13/8)+(-8)*t)+(-772/5)*cos((59/16)+(-7)*t)+(-11115/23)*cos((117/50)+(-6)*t)+(-4310/13)*cos((19/8)+(-5)*t)+(-19665/34)*cos((157/78)+(-3)*t)+(-6949/13)*cos((83/27)+(-1)*t)+(-6805/27)*cos((1/145)+2*t)+(-20828/37)*cos((49/74)+4*t)+(-16299/58)*cos((55/43)+9*t)+(-1449/20)*cos((2/23)+23*t)+(-7448/177)*cos((13/18)+28*t)+(-62/21)*cos((7/16)+31*t));


## Generate node set given x and y coordinate sets
function dist(x, y)
    return [x'; y']
end

## Functions for generating approximate normals using a sorted node1
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
function knnFull(nodes, n)
    kdtree = KDTree(nodes);
    idx, tmp = knn(kdtree, nodes, n, true);
    
    return idx
end

## Functions for setting up RBF envrionment
# Center main node
function center(nodes)
    nodes = nodes .- nodes[:,1];

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


## RBF related functions

# RBF definitons
# 1D
# Gaussian
# ϕ(x1,x2,m) = exp(-abs2(m*(x1-x2)));
# ϕ_x(x1,x2,m) = 2*(m^2)*(x2-x1)*exp(-abs2(m*(x1-x2)));

# PHS
ϕ(x1,x2,m) = abs(x1-x2)^m;
ϕ_x(x1,x2,m) = m*sign(x1-x2)*abs(x1-x2)^(m-1);

# Interpolation weights for 1D
function interpolate(nodes, m)
    n = size(nodes, 2);
    A = zeros(n,n);
    x = nodes[1,:];
    f = nodes[2,:];
    
    # Define A
    for i ∈ 1:n, j ∈ 1:n
        A[i,j] = ϕ(x[i], x[j], m);
    end

    # Solve for weights
    λ = A\f;

    return λ
end

# Interpolation function
function S(t, x, λ, m)
    n = size(x,1);
    nt = size(t,1);
    
    s = zeros(nt);
    for i ∈ 1:n, j ∈ 1:nt
        s[j] += λ[i]*ϕ(t[j], x[i], m);
    end
    
    return s
end

# Interpolation derivative function
function S_x(x, λ, m)
    n = size(x,1);
    
    s = 0;
    for i ∈ 1:n
        s += λ[i]*ϕ_x(0, x[i], m);
    end
    
    return s
end


# Function for tensor normals
function findNormals(nodes, n, m)
    # Find number of nodes
    N = size(nodes, 2);
    
    # Compute approximate normals
    appNorms = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFull(nodes, n);
    
    cent = zeros(2,n);
    rot = cent;
    normals = zeros(2,N);
    λ = zeros(n);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Compute RBF weights
        λ = interpolate(rot, m);

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
    
    return normals
end

## Analysis functions
# Vector error
function vecError(trues, calcs)
    N = size(trues,2);
    normErrs = zeros(N);
    for j ∈ 1:N
        magDiff = norm(trues[:,j]-calcs[:,j]);
        mag = norm(trues[:,j]);
        normErrs[j] = magDiff/mag;
    end
    err = maximum(normErrs);
    return err
end

## Plot functions
# Plot of nodes and vectors
function vectorPlot(points, direction)
    N = size(points, 2);
    vecs = points + 0.2*direction;
    a = scatter(points[1,:],points[2,:],
                aspectratio = :equal,
                legend = false,
                show = false,
                markersize = 2,
                markerstrokealpha = 0,
                markeralpha = .9);
    for i ∈ 1:N
        plot!([points[1,i], vecs[1,i]],
              [points[2,i], vecs[2,i]],
              color = :red,
              legend = false,
              show = false,
              lw = 1)
    end
    return a
end

# Nodes plot
function nodePlot(set1)
    a = scatter(set1[1,:],set1[2,:],
                color = :blue,
                show = false,
                legend = false,
                aspectratio = :equal,
                markersize = 1,
                markerstrokealpha = 0)
    return a
end

# RBF interpolation plot
function interPlot(nodes, λ, m)
    x = nodes[1,:];
    y = nodes[2,:];
    t = range(minimum(x), maximum(x), length = 100);
    s = S(t, x, λ, m);
    
    a = scatter(x,y,
                color = :blue,
                aspectratio = :equal,
                show = false,
                legend = false,
                markersize = 1.5,
                markerstrokealpha = 0)
    plot!(t,s,lw=1)

    return a
end

# Animated process plot
function aniPlot(nodes, n = 10, m = 3, delay = 0.001)
    # Find number of nodes
    N = size(nodes, 2);
    
    # Compute approximate normals
    nmls = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFull(nodes, n);
    
    cent = zeros(2,n);
    rot = cent;
    λ = zeros(n);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(nmls[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Compute RBF weights
        λ = interpolate(rot, m);

        # Plot data
        l = @layout [a b c; d e]
        a = plot(nodePlot(nodes),
                 title = "Initial Node Set",
                 xlims = (-1.1,1.1),
                 ylims = (-1.1,1.1))
        b = plot(nodePlot(nodes[:,idx[i]]),
                 title = "Nearest Neighbors",
                 xlims = (-1.1,1.1),
                 ylims = (-1.1,1.1))
        c = plot(nodePlot(cent),
                 title = "Centered",
                 xlims = (-.25,.25),
                 ylims = (-.25,.25))
        d = plot(nodePlot(rot),
                 title = "Rotated",
                 xlims = (-.25,.25),
                 ylims = (-.25,.25))
        e = plot(interPlot(rot, λ, m),
                 title = "Interpolated",
                 xlims = (-.25,.25),
                 ylims = (-.25,.25))

        f = plot(a,b,c,d,e,
                 layout = l)
        
        display(f)
        sleep(delay)
    end
end

## Main function for finding normals
function comp(N=100, n=10, m1=.1, m2=0.5)
    # Parameterizing our curve
    t = range(0,2*π-2*π/N, length = N);
    
    # Compute true normals
    truNorms = [piNormsX.(t)'; piNormsY.(t)'];
    for i ∈ 1:N
        truNorms[:,i] = truNorms[:,i]/norm(truNorms[:,i]);
    end
    
    # Generate nodes
    nodes = dist(piX.(t), piY.(t))
    appnorms = approxNormals(nodes)
    for i ∈ 1:N
        appnorms[:,i] = appnorms[:,i]/norm(appnorms[:,i]);
    end
    
    normals = findNormals(nodes, n, m1);

    # aniPlot(nodes, n, m1, 0.01)
    
    #a = vectorPlot(nodes, normals);
    # b = vectorPlot(nodes, appnorms);
    #display(a)
    
    a = vecError(truNorms, normals)
    # a = vecError(truNorms, appnorms

    # println(" normals = $a")
    # println("appnorms = $b")
    return a
end

function errs(n,m)
    #neighbors = 3:9;
    nodes = 60:5000:100000;

    a = plot()
    # for n ∈ neighbors
        err = [];
        for N ∈ nodes
            push!(err, log10(comp(N, n, m)));
        end

        a = plot!(log10.(nodes),err,
                  legend = false);
        display(a)
    # end
end

errs(15,5)
