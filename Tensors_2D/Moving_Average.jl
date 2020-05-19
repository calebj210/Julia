## A routine for smoothing our surface vectors using a moving average
## The goal is to increase stability with a small decrease in accuracy.

using EvolvingCurves
using Analysis
using Plots
using LaTeXStrings
include("/home/cajacobs/Documents/Julia/Tensors_2D/Packages/EvolvingCurves.jl/src/Evolving_Backend/Parameterization.jl")


# Function defintion
function movingAvg(nodes, vecs; n = 3)
    N = size(nodes,2)
    
    # Compute nearest neighbor indices
    idx = knnFull(nodes, n)

    # Preallocate space for averged vectors
    avgVecs = zeros(2,N)
    
    # Begin moving average using KNN to nodes to average
    for i∈1:N
        avgVecs[:,i] = sum(vecs[:,idx[i]], dims = 2)
    end

    return avgVecs/n
end


### Test space
## Circle parameterization
# Unit circle
x(t,r) = r*cos.(t)
y(t,r) = r*sin.(t)

## Growing Circle
function grow(N=20, n=5, m=3, o=0; r=.2, TF = 0.02, Δt = 10^(-4))
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N)
    
    # Creating our node distribution
    nodes = [x(t,r)';y(t,r)']

    # Compute time partition
    time = Δt:Δt:TF

    # Initiate error variables
    i = 1
    radiusErr = zeros(size(time,1));
    compR = zeros(N);
    
    # Begin evolving our curve
    for T∈time
        # Compute normal direction at each node
        norms = findNormals(nodes, n, m, o);
        if norms[:,1]⋅[1,0] < 0
            norms *= -1
        end

        # Compute curvature at each node
        κ = computeκ(nodes, n, m, o)
        κ = [κ';κ']

        # Compute curvature scaled normal at each node
        norms = κ.*norms

        # Compute moving average of nodes
        norms = movingAvg(nodes, norms, n = 11)

        # Compute the position of each node at the next tine step
        nodes += norms*Δt

        # Plot evolving surface
        # plotA = scatter(nodes[1,:],nodes[2,:],
        #                 color = :blue,
        #                 legend = false,
        #                 aspectratio = :equal,
        #                 markersize = 1,
        #                 markerstrokealpha = 0,
        #                 xlims = (-0.3,0.3),
        #                 ylims = (-0.3,0.3));
        # display(plotA)
        # # sleep(0.01)

        # Compute largest error at each time step
        for j ∈ 1:N
            compR[j] = norm(nodes[:,j]);
        end
        truR = fill((r^2+2T)^(0.5),N);
        radiusErr[i] = norm(scalError(truR,compR)[2],Inf); 
        i += 1;
    end

    # Plot the largest error over time
    plotB = plot(time,radiusErr,
                 title = "Relative Error for Growing Circle",
                 xlabel = "Time (s)",
                 ylabel = L"||\mathrm{Rel}\;\mathrm{Err}||_\infty",
                 yscale = :log10,
                 xlims = (0,TF),
                 xticks = [0:TF/10:TF...])
    display(plotB)
end

grow(60, 5, 3, 3, r = 0.2, TF = 0.02, Δt = 10^(-4))
