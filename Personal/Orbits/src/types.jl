# src/types.jl

"""
    Body

Represents a physical body in the 2D gravitational simulation.

# Fields
- `name::String`: Name of the body (used for labels).
- `mass::Float64`: Mass of the body.
- `x::Float64`: Initial x-coordinate of position.
- `y::Float64`: Initial y-coordinate of position.
- `vx::Float64`: Initial x-coordinate of velocity.
- `vy::Float64`: Initial y-coordinate of velocity.
- `color::Any`: Color used for visualization (e.g., Symbol like `:blue`, String, or RGBA/RGB from Makie).
- `size::Float64`: Visual marker size of the body.
"""
struct Body
    name::String
    mass::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
    color::Any
    size::Float64
end

"""
    Body(; name="Body", mass=1.0, pos=(0.0, 0.0), vel=(0.0, 0.0), color=:black, size=15.0)

Convenient keyword constructor for `Body`.
"""
function Body(;
    name::String="Body",
    mass::Real=1.0,
    pos::Tuple{Real, Real}=(0.0, 0.0),
    vel::Tuple{Real, Real}=(0.0, 0.0),
    color=Symbol(:black),
    size::Real=15.0
)
    return Body(name, Float64(mass), Float64(pos[1]), Float64(pos[2]), Float64(vel[1]), Float64(vel[2]), color, Float64(size))
end
