### Computing shrinking and growing circles
using Plots
using EvolvingCurves
using Analysis
using LinearAlgebra
using LaTeXStrings
include("Node_Adapt.jl")
gr()

## Unit circle
x(r,t) = r*cos.(t);
y(r,t) = r*sin.(t);

## Small Circle
x2(t) = 0.2*cos.(t) .+ 0.4;
y2(t) = 0.2*sin.(t) .+ 0.4;

## True normals
truX(t) = cos.(t);
truY(t) = sin.(t);

## Curvature testing
function main1(N=20, n=5, m=3, o=0)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N);

    # Creating our node distribution
    nodes = [x(t)';y(t)'];

    κ1 = computeK(nodes, n, m, o)[1] - 1;
    display(κ1)

    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = 2N);

    # Creating our node distribution
    nodes = [x2(t)';y2(t)'];

    κ2 = computeK(nodes, n, m, o)[1] - 1;
    display(κ2)
    
    display(log2(abs(κ1/κ2)));
end

## 1 Circles
function main2(N=20, n=5, m=3, o=0, Δt = 0.01, TF=0.02, r=0.2)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N);

    # Creating our node distribution
    nodes = [x(r,t)';y(r,t)'];
    
    # Evolve curve
    time = Δt:Δt:TF;
    radiusErr = zeros(size(time,1));
    compR = zeros(N);
    i = 1;
    k = 1;
    T = 0;
    Nodes = nodes[:];
    for T ∈ time
        if k > 10
            Nodes = Nodes[1:2N];
            nodeAdapt!(Nodes, n, m, o);
            nodes = reshape(Nodes,2,:);
            k=1;
            display(size(Nodes,1)/2)
        end
        k += 1;
        
        # Compute normal direction
        norms = findNormals(nodes, n, m, o);
        if norms[:,1]⋅[1,0] < 0
            norms *= -1;
        end
        # Compute curvature
        κ = computeκ(nodes, n, m, o);
        K = [κ';κ']

        # display(κ[1])
        
        nodes += K.*norms*Δt;
        
        # plotA = scatter(nodes[1,:],nodes[2,:],
        #                 color = :blue,
        #                 show = false,
        #                 legend = false,
        #                 aspectratio = :equal,
        #                 markersize = 1,
        #                 markerstrokealpha = 0,
        #                 xlims = (-0.3,0.3),
        #                 ylims = (-0.3,0.3));
        # # display(T)
        # display(plotA)
        # sleep(0.01)
        
        for j ∈ 1:N
            compR[j] = norm(nodes[:,j]);
        end
        truR = fill((r^2+2T)^(0.5),N);
        radiusErr[i] = norm(scalError(truR,compR)[2],Inf); 
        i += 1;
    end
    # plotA = scatter(nodes[1,:],nodes[2,:],
    #                 color = :blue,
    #                 show = false,
    #                 legend = false,
    #                 aspectratio = :equal,
    #                 markersize = 2,
    #                 markerstrokealpha = 0,
    #                 xlims = (-2,2),
    #                 ylims = (-2,2));
    plotB = plot(time,radiusErr,
                 title = "Relative Error for Growing Circle",
                 xlabel = "Time (s)",
                 ylabel = L"||\mathrm{Rel}\;\mathrm{Err}||_\infty",
                 yscale = :log10,
                 xlims = (0,TF),
                 xticks = [0:TF/10:TF...])
    display(plotB)

    
    # return radiusErr
end

## 2 Circles
function main3(N=20, n=5, m=3, o=0, Δt = 0.01)
    NN = round(Int,N/2);
    
    # Coordinate discretization
    t = range(0,2*π-2*π/NN, length = NN);

    # Creating our node distribution
    nodes = zeros(2,N)
    nodes[:,1:NN] = [x2(t)';y2(t)'];
    nodes[:,NN+1:end] = [x2(t)';y2(t)'] .+ 0.2;
    
    # Evolve curve                
    T = Δt;
    while T < 0.01
        # Compute normal direction
        norms = findNormals(nodes, n, m, o)

        # Compute curvature
        κ = computeκ(nodes, n, m, o);
        K = [κ';κ']
        
        nodes += K.*norms*Δt;

        T += Δt;
        
        plotA = scatter(nodes[1,:],nodes[2,:],
                        color = :blue,
                        show = false,
                        legend = false,
                        aspectratio = :equal,
                        markersize = 2,
                        markerstrokealpha = 0,
                        xlims = (0.2,0.8),
                        ylims = (0.2,0.8));
        display(T)
        display(plotA)
        # sleep(0.01)
    end
end

## Complete Error Plots
# Error with varying the total number of nodes
function fig_a()
    # Number of nodes to analyze
    nodes = [50,100,150,200,300]

    # Time interval
    t = 10^(-6):10^(-6):0.02;

    # Prelocate plot
    A = plot();
    
    # Produce plots
    for N ∈ nodes
        errs = main2(N,12,5,5,10^(-6),0.02);
        A = plot!(t,errs,
                  title = "Relative Error for Shrinking Circle",
                  xlabel = "Time (s)",
                  ylabel = L"||\mathrm{Rel}\;\mathrm{Err}\;||_\infty",
                  yscale = :log10,
                  xlims = (0,0.02),
                  label = string(N," Nodes"),
                  legend = :topleft,
                  dpi=300);
        display(A)
    end

    # Export plot to fig_a.png
    savefig(A,"fig_a.png")
end 

# Error with varying the number of neighbors
function fig_b()
    # Number of neigbors
    neigbs = [6,8,10,12,14]

    # Time interval
    t = 10^(-6):10^(-6):0.02;

    # Prelocate plot
    B = plot();
    
    # Produce plots
    for N ∈ neigbs
        errs = main2(100,N,5,floor(Int,N/2)-1,10^(-6),0.02);
        B = plot!(t,errs,
                  title = "Relative Error for Shrinking Circle",
                  xlabel = "Time (s)",
                  ylabel = L"||\mathrm{Rel}\;\mathrm{Err}\;||_\infty",
                  yscale = :log10,
                  xlims = (0,0.02),
                  label = string(N," Neighbors"),
                  legend = :topleft,
                  dpi=300);
        display(B)
    end

    # Export plot to fig_a.png
    savefig(B,"fig_b.png")
end


main2(50, 3, 3, 2, 10^(-4), 0.02,0.2)
