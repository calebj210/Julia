### Laplace Beltrami test on an ellipse
using Plots
using EvolvingCurves
using Analysis
using LaTeXStrings
gr()

## Ellipse parameterization
x(a,t) = a*cos.(t)
y(b,t) = b*sin.(t)

## Surface function
f(x,y,Pp,Pq) = (x.^Pp)*(y.^Pq)

## True LBO
lap_R2_X(x,y,Pp,Pq) = (Pp < 2) ? (0) : ((-1+Pp)*Pp*x.^(Pp-2).*y.^Pq)
lap_R2_Y(x,y,Pp,Pq) = (Pq < 2) ? (0) : ((-1+Pq)*Pq*x.^Pp.*y.^(Pq-2))
trueLBO(x,y,Pp,Pq,a,b) =
    lap_R2_X(x,y,Pp,Pq) + lap_R2_Y(x,y,Pp,Pq) -
    x.^Pp.*y.^Pq.*(a^4*b^4*(b^2*Pp+a^2*Pq))./(b^4*x.^2+a^4*y.^2).^2-(b^4*(-1+Pp)*Pp+2*a^2*b^2*Pp*Pq+a^4*(Pq-1)*Pq)*x.^Pp.*y.^Pq./(b^4*x.^2+a^4*y.^2)


## Test function for computing errors
function computeErrs(NN=[20]; n=5, m=3, o=0, a=2, b=1, Pp=0, Pq=2)
    # Initialize error storage
    errs = zeros(size(NN,1))
    
    for j = 1:size(NN,1)
        # Get the number of nodes
        N = NN[j]
        
        # Coordinate discretization
        t = range(0,2*π-2*π/N, length = N)
        
        # Create our node distribution
        nodes = [x.(a,t)';y.(b,t)']
        
        # Compute function values at each node
        F = zeros(N)
        for i∈1:N
            F[i] = f(nodes[1,i],nodes[2,i],Pp,Pq)
        end
        
        # Discretize LBO
        D = constructLBO(nodes, n, m, o)
        
        # Compute the LBO
        LBOF = D*F
        
        # Compute true solution
        trueLBOF = zeros(N)
        for i∈1:N
            trueLBOF[i] = trueLBO(nodes[1,i],nodes[2,i],Pp,Pq,a,b)
        end
        
        # Compute errors
        errs[j] = scalError(trueLBOF,LBOF)[1]
    end

    return errs
end

# Convergence plot function
function convPlot(;NN = [20,40,80,160,320,640,1280], O=[2])
    # Initialize the convergence plot
    plotA = plot()

    # Loop through polynomial degrees
    for o∈O
        # Compute errors at specified polynomial degree
        errs = computeErrs(NN, n=2(o+1), m=5, o=o, Pp=0, Pq=2, a=2, b=1)

        plot!(NN,errs,
              title = "Convergence of LBO on an Ellipse",
              xlabel = "# Nodes",
              ylabel = L"||\mathrm{Rel}\;\mathrm{Err}\;||_\infty",
              xscale = :log10,
              yscale = :log10,
              label = string("Degree ", o, " Polynomial"),
              legend = :topright,
              dpi = 300)
    end

    # Export convergence to plot a png file
    savefig(plotA,"Convergence.png")

    # Display converence plot
    display(plotA)
end

## Run test
# computeErrs([20,40,80,160,320,640,1280], n=16, m=5, o=7)

## Genergerate convergence plot
convPlot(NN=[20:1000:20000...], O=[2,3,4,5,6,7])
