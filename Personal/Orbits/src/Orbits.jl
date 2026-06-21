module Orbits

# Include component files
include("types.jl")
include("simulation.jl")
include("visualization.jl")

# Export types
export Body

# Export simulation functions
export simulate_orbits

# Export visualization functions
export visualize_orbits, save_orbit_animation, simulate_and_visualize

"""
    simulate_and_visualize(bodies::Vector{Body}, tspan::Tuple{Real, Real}; G=1.0, softening=0.0, solver=Tsit5(), reltol=1e-8, abstol=1e-8, trail_duration=nothing, title="2D Gravitational Orbits", xlims=nothing, ylims=nothing, kwargs...)

Simulates the orbits of `bodies` over `tspan` and opens the interactive visualization dashboard.

See `simulate_orbits` and `visualize_orbits` for arguments and details.
"""
function simulate_and_visualize(
    bodies::Vector{Body},
    tspan::Tuple{Real, Real};
    G::Real=1.0,
    softening::Real=0.0,
    solver=Tsit5(),
    reltol=1e-8,
    abstol=1e-8,
    trail_duration::Union{Nothing, Real}=nothing,
    title::String="2D Gravitational Orbits",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    kwargs...
)
    sol = simulate_orbits(
        bodies,
        tspan;
        G = G,
        softening = softening,
        solver = solver,
        reltol = reltol,
        abstol = abstol,
        kwargs...
    )
    
    return visualize_orbits(
        sol,
        bodies;
        trail_duration = trail_duration,
        title = title,
        xlims = xlims,
        ylims = ylims
    )
end

end # module Orbits
