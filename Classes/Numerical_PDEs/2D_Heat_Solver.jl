### 2D Heat Equation FD Solver
# Loading packages
using LinearAlgebra
using Plots
plotly

# Node matrix constructor
function genA(M,N,BC)
    U = zeros(M+2,N+2);
    
    U[1,:] .= BC;
    U[M+2,:] .= BC;
    U[:,1] .= BC;
    U[:,N+2] .= BC;
    
    return U
end

# Time stepper
function tStep(U0,dt,X,Y)
    dx = (X[:b] - X[:a])/(X[:M]+1);
    dy = (Y[:b] - Y[:a])/(Y[:N]+1);

    # Preallocate space for next step
    U = deepcopy(U0);

    # Begin updating U
    for i ∈ 2:X[:M]+1
        for j ∈ 2:Y[:N]+1
            U[i,j] =
                U0[i,j] + dt*((U0[i,j+1] - 2U0[i,j] + U0[i,j-1])/(dx^2)
                              +(U0[i+1,j] - 2U0[i,j] + U0[i-1,j])/(dy^2));
        end
    end

    return U
end

# Complete solver
function comp(ax,bx,ay,by,M,N,dt,BC)
    U = genA(M,N,BC);
    U[5:10,5:10] .= 1000;
    U[20:25,20:25] .= 1000;

    X = Dict([(:a,ax),(:b,bx),(:M,M)]);
    Y = Dict([(:a,ay),(:b,by),(:N,N)]);
    
    for i ∈ 1:1000
        U = tStep(U,dt,X,Y);
        
        display(heatmap(U,
                        clims = (0,1000)))
        sleep(0.01)
    end
end

comp(0,1,0,1,30,30,0.00001,500)
