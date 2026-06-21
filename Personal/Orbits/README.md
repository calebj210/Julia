# Orbits.jl

`Orbits.jl` is a lightweight, high-performance 2D gravitational N-body orbit simulator and interactive visualization tool written in Julia. It integrates standard differential equation solvers with interactive graphics to animate complex orbital systems.

It relies on **OrdinaryDiffEq.jl** (part of the `DifferentialEquations.jl` ecosystem) for adaptive time stepping and ODE evolution, and **Makie.jl** for premium interactive dashboards and video rendering.

---

## 🌌 Mathematical Model

The package simulates $N$ point-mass bodies moving in a 2-dimensional plane. Under Newton's law of universal gravitation, the equations of motion for each body $i \in \{1, \dots, N\}$ with mass $m_i$ and coordinates $\mathbf{q}_i = (x_i, y_i)$ are:

$$\frac{d\mathbf{q}_i}{dt} = \mathbf{v}_i$$

$$\frac{d\mathbf{v}_i}{dt} = \sum_{j \neq i} G m_j \frac{\mathbf{q}_j - \mathbf{q}_i}{\left( \|\mathbf{q}_j - \mathbf{q}_i\|^2 + \epsilon^2 \right)^{3/2}}$$

Where:
- $G$ is the gravitational constant (default: `1.0`).
- $\epsilon \ge 0$ is an optional **softening parameter** used to prevent numerical singularities (division by zero) during close encounters or collisions.

---

## 🚀 Installation

Ensure you have Julia (v1.9+) installed. Since this package is located in your local Orbits folder, you can activate and run it as follows:

```julia
using Pkg
# Activate this local package environment
Pkg.activate(".")
Pkg.instantiate()
```

---

## 🛠️ API Reference

### 1. `Body`
Represents an orbiting body.
```julia
Body(; name="Body", mass=1.0, pos=(x, y), vel=(vx, vy), color=:black, size=15.0)
```
- Defined in [src/types.jl](file:///home/merlin/Documents/Julia/Personal/Orbits/src/types.jl).
- `color` can be a `Symbol` (e.g. `:blue`), a string (e.g. `"red"`), or RGB/RGBA values.
- `size` controls the radius of the marker in the plot.

### 2. `simulate_orbits`
Integrates the ODE equations of motion.
```julia
sol = simulate_orbits(bodies, tspan; G=1.0, softening=0.0, solver=Tsit5(), reltol=1e-8, abstol=1e-8, kwargs...)
```
- Defined in [src/simulation.jl](file:///home/merlin/Documents/Julia/Personal/Orbits/src/simulation.jl).
- Returns a standard `ODESolution` from `OrdinaryDiffEq`.
- Because the solvers use **adaptive time stepping**, the solution object contains an $O(1)$ continuous interpolant. We exploit this during animation to sample coordinates on a perfectly uniform grid, keeping playback speeds constant.

### 3. `visualize_orbits`
Launches the interactive dashboard GUI.
```julia
fig = visualize_orbits(sol, bodies; trail_duration=nothing, title="2D Orbits", xlims=nothing, ylims=nothing)
```
- Defined in [src/visualization.jl](file:///home/merlin/Documents/Julia/Personal/Orbits/src/visualization.jl).
- Returns a Makie `Figure` object.
- Automatically computes stable axis limits bounding the entire trajectory.
- Features: Play/Pause button, Time scrubber slider, Playback Speed scale slider, Reset button, and Fading/Finitude Trail toggle.

### 4. `save_orbit_animation`
Renders and saves the simulation to a file.
```julia
save_orbit_animation(sol, bodies, "filename.gif"; fps=30, duration=10.0, trail_duration=nothing, xlims=nothing, ylims=nothing)
```
- Defined in [src/visualization.jl](file:///home/merlin/Documents/Julia/Personal/Orbits/src/visualization.jl).
- Supports `.gif` and `.mp4` formats.

### 5. `simulate_and_visualize`
A wrapper that simulates and immediately launches the visual dashboard.
```julia
simulate_and_visualize(bodies, tspan; G=1.0, solver=Tsit5(), ...)
```
- Defined in [src/Orbits.jl](file:///home/merlin/Documents/Julia/Personal/Orbits/src/Orbits.jl).

---

## 🎨 Interactive GUI Features

When you call `visualize_orbits(sol, bodies)` or `simulate_and_visualize(...)`, an interactive GUI opens:

```
  +-------------------------------------------------------------+
  |              2D Gravitational Orbits Simulation             |
  |                                                             |
  |      [Body A]       *---..                                  |
  |      (scatter)      \      `---.                            |
  |                      \          `* [Body B]                 |
  |                       \          /                          |
  |                        *________/                           |
  |                     [Body C]                                |
  |                                                             |
  +-------------------------------------------------------------+
  |  [==== Slider: Time ===============================]  Time: 4.2 |
  |  [Play]  [Reset]   Playback Speed: [=====] 1.00x            |
  +-------------------------------------------------------------+
```

1. **Play / Pause**: Starts or pauses the orbital evolution. Runs asynchronously in a background loop without blocking your Julia REPL.
2. **Time Scrubber Slider**: Click or drag to jump to any point in the simulation.
3. **Playback Speed Slider**: Speed up or slow down playback scale dynamically.
4. **Reset**: Instantly pauses the playback and returns all bodies to $t = 0$.
5. **Dynamic Fading Trail**: The `trail_duration` parameter limits how much orbital history is visible, drawing a sliding tail behind the body.

---

## 📖 Quick Examples

### 1. Simple Keplerian Orbit (Sun + Earth)

```julia
using Orbits

# Define a massive stationary Sun and a light orbiting Earth
sun = Body(name="Sun", mass=1.0, pos=(0.0, 0.0), vel=(0.0, 0.0), color=:yellow, size=25.0)
earth = Body(name="Earth", mass=1e-6, pos=(1.0, 0.0), vel=(0.0, 1.0), color=:blue, size=12.0)

# Simulate for one year (T = 2π)
fig = simulate_and_visualize([sun, earth], (0.0, 2*pi); G=1.0, trail_duration=1.5)

# In GLMakie, returning the figure will display the interactive window:
display(fig)
```

### 2. Choreographed Figure-8 Orbit (Three Equal Masses)

See the full script in [examples/figure_eight.jl](file:///home/merlin/Documents/Julia/Personal/Orbits/examples/figure_eight.jl).

```julia
using Orbits
using OrdinaryDiffEq

# Initial positions and velocities for figure-8 choreographic orbit
x1, y1 = -0.97000436, 0.24308753
vx1, vy1 = 0.46620531, 0.43236573

b1 = Body(name="A", mass=1.0, pos=(x1, y1), vel=(vx1, vy1), color=:cyan, size=15)
b2 = Body(name="B", mass=1.0, pos=(-x1, -y1), vel=(vx1, vy1), color=:magenta, size=15)
b3 = Body(name="C", mass=1.0, pos=(0.0, 0.0), vel=(-2*vx1, -2*vy1), color=:orange, size=15)

# Simulate for 3 orbital periods
t_period = 6.3259
fig = simulate_and_visualize([b1, b2, b3], (0.0, 3 * t_period), G=1.0, trail_duration=2.0)

display(fig)
```

### 3. Save to a GIF headlessly

```julia
using Orbits

# Setup bodies
sun = Body(name="Sun", mass=1.0, pos=(0.0, 0.0), vel=(0.0, 0.0), color=:orange, size=25)
earth = Body(name="Earth", mass=1e-5, pos=(1.0, 0.0), vel=(0.0, 1.0), color=:blue, size=12)

# Solve ODE
sol = simulate_orbits([sun, earth], (0.0, 10.0); G=1.0)

# Save animation without opening a window (uses precalculated bounding box)
save_orbit_animation(sol, [sun, earth], "earth_orbit.gif", fps=30, duration=8.0, trail_duration=2.0)
```
