### Circle in a Vortex

## Including used packages
import Analysis
using Plots
gr()

## Problem defintitions
# Unit circle
x(θ) = 0.15cos.(θ) .+ 0.5;
y(θ) = 0.15sin.(θ) .+ 0.75;

# Velocity field
xVel(x,y,t) = -sin.(π*x).^2*sin.(2*π*y)*cos.(π*t/4);
yVel(x,y,t) = sin.(2*π*x)*sin.(π*y).^2*cos.(π*t/4);

## Solver and Evolver
function main(N=20, n=5, m=3, o=0, Δt=0.1)
    # Coordinate discretization
    θ = range(0,2*π-2*π/N, length = N);

    # Creating our node distribution
    nodes = [x(θ)';y(θ)'];

    t = 0;
    while t < 4-Δt
        for i ∈ 1:N
            nodes[:,i] += [xVel(nodes[1,i],nodes[2,i],t)'; yVel(nodes[1,i],nodes[2,i],t)']*Δt;
        end
        
        t += Δt;
        
        plotA = scatter(nodes[1,:],nodes[2,:],
                        color = :blue,
                        show = false,
                        legend = false,
                        aspectratio = :equal,
                        markersize = 1.5,
                        markerstrokealpha = 0,
                        xlims = (0,1),
                        ylims = (0,1));
        display(t)
        display(plotA)
        sleep(Δt)
    end
end

main(100,5,3,2,0.01)
