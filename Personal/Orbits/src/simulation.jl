# src/simulation.jl

using OrdinaryDiffEq

"""
    gravity_ode!(du, u, p, t)

The ODE system representing Newtonian gravity in 2D for N bodies.
`u` is a flat vector of length 4N structured as:
[x1, y1, vx1, vy1, x2, y2, vx2, vy2, ...]
`p` is a tuple containing `(masses, G, softening)`.
"""
function gravity_ode!(du, u, p, t)
    masses, G, softening = p
    N = length(masses)
    
    # 1. Update positions (dx/dt = vx, dy/dt = vy)
    for i in 1:N
        idx_pos = 4*(i-1) + 1
        idx_vel = 4*(i-1) + 3
        
        du[idx_pos]   = u[idx_vel]     # dx/dt = vx
        du[idx_pos+1] = u[idx_vel+1]   # dy/dt = vy
    end
    
    # 2. Update velocities (dvx/dt = ax, dvy/dt = ay)
    for i in 1:N
        idx_pos_i = 4*(i-1) + 1
        idx_vel_i = 4*(i-1) + 3
        
        xi = u[idx_pos_i]
        yi = u[idx_pos_i+1]
        
        ax = 0.0
        ay = 0.0
        
        for j in 1:N
            i == j && continue
            idx_pos_j = 4*(j-1) + 1
            xj = u[idx_pos_j]
            yj = u[idx_pos_j+1]
            
            dx = xj - xi
            dy = yj - yi
            r2 = dx^2 + dy^2 + softening^2
            r = sqrt(r2)
            
            # Acceleration contribution: G * m_j * dx / r^3
            if r > 1e-12  # Avoid division by zero
                factor = G * masses[j] / (r * r2)
                ax += factor * dx
                ay += factor * dy
            end
        end
        
        du[idx_vel_i]   = ax
        du[idx_vel_i+1] = ay
    end
    
    return nothing
end

"""
    simulate_orbits(bodies::Vector{Body}, tspan::Tuple{Real, Real}; G=1.0, softening=0.0, solver=Tsit5(), reltol=1e-8, abstol=1e-8, kwargs...)

Simulates the motion of `bodies` over `tspan` using `OrdinaryDiffEq`.

# Arguments
- `bodies`: A vector of `Body` instances defining the initial state.
- `tspan`: A tuple `(t_start, t_end)` for the integration time interval.
- `G`: Gravitational constant (default: 1.0).
- `softening`: A small softening parameter to avoid singularities on collisions (default: 0.0).
- `solver`: The ODE solver from OrdinaryDiffEq (default: `Tsit5()`).
- `reltol`, `abstol`: Relative and absolute solver tolerances (default: 1e-8).
- `kwargs...`: Additional keyword arguments passed to `solve`.

# Returns
- A `ODESolution` object which can be evaluated at any time `t` in `tspan`.
"""
function simulate_orbits(
    bodies::Vector{Body},
    tspan::Tuple{Real, Real};
    G::Real=1.0,
    softening::Real=0.0,
    solver=Tsit5(),
    reltol=1e-8,
    abstol=1e-8,
    kwargs...
)
    N = length(bodies)
    
    # Pack initial conditions into a flat vector
    u0 = zeros(4 * N)
    for i in 1:N
        b = bodies[i]
        idx_pos = 4*(i-1) + 1
        idx_vel = 4*(i-1) + 3
        
        u0[idx_pos]     = b.x
        u0[idx_pos+1]   = b.y
        u0[idx_vel]     = b.vx
        u0[idx_vel+1]   = b.vy
    end
    
    # Parameters: masses, G, softening
    masses = [b.mass for b in bodies]
    p = (masses, Float64(G), Float64(softening))
    
    tspan_float = (Float64(tspan[1]), Float64(tspan[2]))
    
    prob = ODEProblem(gravity_ode!, u0, tspan_float, p)
    sol = solve(prob, solver; reltol=reltol, abstol=abstol, kwargs...)
    
    return sol
end
