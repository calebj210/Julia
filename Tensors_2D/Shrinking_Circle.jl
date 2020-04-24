### Shrinking Circle

## Including used packages
import EvolvingCurves.findNormals
import Analysis
using Plots
gr()

## Surface Definition
# Unit circle
x(t) = cos.(t);
y(t) = sin.(t);

# True solution
truX(θ,t) = (1-t)cos.(θ);
truY(θ,t) = (1-t)sin.(θ);

## Solver and Evolver
function main(N=20, n=5, m=3, o=0, Δt=0.1)
    # Coordinate discretization
    t = range(0,2*π-2*π/N, length = N);

    # Creating our node distribution
    nodes = [x(t)';y(t)'];

    T = Δt;
    while T < 1
        norms = -findNormals(nodes, n, m, o)
        nodes += norms*Δt;

        T += Δt;
        
        plotA = scatter(nodes[1,:],nodes[2,:],
                        color = :blue,
                        show = false,
                        legend = false,
                        aspectratio = :equal,
                        markersize = 2,
                        markerstrokealpha = 0,
                        xlims = (-1,1),
                        ylims = (-1,1));
        display(T)
        display(plotA)
        sleep(Δt)
    end
end

main(20,5,3,2,0.1)
