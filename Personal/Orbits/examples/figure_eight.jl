# examples/figure_eight.jl

using Orbits
using OrdinaryDiffEq

# --- Define the Figure-8 Initial Conditions ---
# The Chenciner-Montgomery figure-eight orbit is a choreographic solution
# to the three-body problem where three equal masses (m = 1.0) move along
# the same path in a plane. G = 1.0.

x1, y1 = -0.97000436, 0.24308753
vx1, vy1 = 0.46620531, 0.43236573

# By symmetry:
# Body 1
b1 = Body(
    name = "Body A",
    mass = 1.0,
    pos = (x1, y1),
    vel = (vx1, vy1),
    color = :deepskyblue,
    size = 14
)

# Body 2
b2 = Body(
    name = "Body B",
    mass = 1.0,
    pos = (-x1, -y1),
    vel = (vx1, vy1),
    color = :crimson,
    size = 14
)

# Body 3 (at origin, with balancing momentum)
b3 = Body(
    name = "Body C",
    mass = 1.0,
    pos = (0.0, 0.0),
    vel = (-2 * vx1, -2 * vy1),
    color = :gold,
    size = 14
)

bodies = [b1, b2, b3]

# Time span: one orbital period is roughly T = 6.3259
period = 6.32591398
tspan = (0.0, 2.5 * period) # Simulate 2.5 periods

println("Running Figure-8 orbit simulation...")
# Use Tsit5 with tight tolerances to preserve the orbit
sol = simulate_orbits(bodies, tspan; G = 1.0, reltol = 1e-10, abstol = 1e-10)

output_file = "figure_eight.gif"
println("Saving animation to $output_file...")
# We use a trail duration of 1.5 time units to show a beautiful fading tail
save_orbit_animation(
    sol,
    bodies,
    output_file;
    fps = 30,
    duration = 12.0, # 12 seconds video duration
    trail_duration = 1.5,
    title = "Periodic 3-Body Figure-8 Orbit",
    xlims = (-1.2, 1.2),
    ylims = (-0.6, 0.6)
)
println("Done! Animation saved successfully.")
