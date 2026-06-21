# src/visualization.jl

using Makie, GLMakie

"""
    get_position_at_t(sol, t::Real, i::Int)

Extracts the position of the `i`-th body at simulation time `t` as a `Point2f`.
"""
function get_position_at_t(sol, t::Real, i::Int)
    state = sol(t)
    return Point2f(state[4*(i-1) + 1], state[4*(i-1) + 2])
end

"""
    get_trail_at_t(sol, t::Real, i::Int, trail_duration::Union{Nothing, Real})

Generates a vector of `Point2f` representing the trail of the `i`-th body up to time `t`.
If `trail_duration` is `nothing`, the trail starts at the beginning of the simulation.
Otherwise, it starts at `t - trail_duration`.
"""
function get_trail_at_t(sol, t::Real, i::Int, trail_duration::Union{Nothing, Real})
    t_start = isnothing(trail_duration) ? sol.t[1] : max(sol.t[1], t - trail_duration)
    
    # If the time span is extremely short, return the current position repeated to make a valid trail line
    if t - t_start < 1e-6
        p = get_position_at_t(sol, t, i)
        return [p, p]
    end
    
    ts = range(t_start, t, length=1000)
    points = Vector{Point2f}(undef, length(ts))
    for (k, tk) in enumerate(ts)
        state = sol(tk)
        points[k] = Point2f(state[4*(i-1) + 1], state[4*(i-1) + 2])
    end
    return points
end

"""
    visualize_orbits(sol, bodies::Vector{Body}; trail_duration=nothing, title="2D Gravitational Orbits", xlims=nothing, ylims=nothing)

Creates an interactive Makie GUI dashboard for exploring the simulated orbits.

# Arguments
- `sol`: The ODESolution returned by `simulate_orbits`.
- `bodies`: The vector of `Body` instances.
- `trail_duration`: How long the orbital tail should be in simulation time units. If `nothing`, shows the full history.
- `title`: The window/axis title.
- `xlims`, `ylims`: Fixed axis limit tuples `(min, max)` to override automatic calculations.

# Returns
- A `Figure` object. Displaying or returning this figure opens the interactive GUI window.
"""
function visualize_orbits(
    sol,
    bodies::Vector{Body};
    trail_duration::Union{Nothing, Real}=nothing,
    title::String="2D Gravitational Orbits",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing
)
    N = length(bodies)
    
    # 1. Pre-calculate default axis limits based on the full simulation trajectory
    xs = Float64[]
    ys = Float64[]
    # Sample 300 points along the path to determine spatial extent
    for t in range(sol.t[1], sol.t[end], length=300)
        state = sol(t)
        for i in 1:N
            push!(xs, state[4*(i-1) + 1])
            push!(ys, state[4*(i-1) + 2])
        end
    end
    
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)
    dx = xmax - xmin
    dy = ymax - ymin
    
    # Handle zero range case (e.g. single stationary body)
    if dx ≈ 0.0; dx = 1.0; end
    if dy ≈ 0.0; dy = 1.0; end
    
    padding = 0.1
    default_xlims = (xmin - padding*dx, xmax + padding*dx)
    default_ylims = (ymin - padding*dy, ymax + padding*dy)
    
    # 2. Setup Figure and Main Axis
    fig = Figure(size = (1000, 850))
    ax = Axis(
        fig[1, 1],
        aspect = DataAspect(),
        title = title,
        xlabel = "x",
        ylabel = "y"
    )
    
    # Apply axis limits
    if !isnothing(xlims)
        xlims!(ax, xlims...)
    else
        xlims!(ax, default_xlims...)
    end
    
    if !isnothing(ylims)
        ylims!(ax, ylims...)
    else
        ylims!(ax, default_ylims...)
    end
    
    # 3. Create Observables for reactive plotting
    t_obs = Observable(sol.t[1])
    
    pos_observables = [
        lift(t -> get_position_at_t(sol, t, i), t_obs) for i in 1:N
    ]
    
    trail_observables = [
        lift(t -> get_trail_at_t(sol, t, i, trail_duration), t_obs) for i in 1:N
    ]
    
    # 4. Plot Trails and Bodies
    for i in 1:N
        color_val = bodies[i].color
        # Convert simple symbol color to tuple color to apply transparency
        trail_color = color_val isa Symbol ? (color_val, 0.6) : color_val
        
        # Plot orbital trail
        lines!(
            ax,
            trail_observables[i],
            color = trail_color,
            linewidth = 2
        )
        
        # Plot body
        scatter!(
            ax,
            pos_observables[i],
            color = color_val,
            markersize = bodies[i].size,
            marker = :circle,
            strokecolor = :white,
            strokewidth = 1.0
        )
        
        # Plot body name label
        if !isempty(bodies[i].name)
            text!(
                ax,
                pos_observables[i],
                text = bodies[i].name,
                align = (:center, :bottom),
                offset = (0, 8),
                fontsize = 12,
                color = :black
            )
        end
    end
    
    # 5. Interactive GUI Controls Panel
    controls = fig[2, 1] = GridLayout()
    
    # Time scrubber slider
    t_slider = Slider(controls[1, 1], range = range(sol.t[1], sol.t[end], length=1000), startvalue = sol.t[1])
    t_label = Label(controls[1, 2], @lift(string("Time: ", round($(t_slider.value), digits=2))), width = 120, halign = :left)
    
    # Connect the slider's value back to our plot's time observable
    connect!(t_obs, t_slider.value)
    
    # Control buttons and playback speed slider
    btn_layout = controls[2, 1] = GridLayout()
    
    play_btn = Button(btn_layout[1, 1], label = "Play", width = 80)
    reset_btn = Button(btn_layout[1, 2], label = "Reset", width = 80)
    
    Label(btn_layout[1, 3], "Playback Speed:")
    speed_slider = Slider(btn_layout[1, 4], range = 0.05:0.05:5.0, startvalue = 1.0, width = 150)
    speed_label = Label(btn_layout[1, 5], @lift(string(round($(speed_slider.value), digits=2), "x")), width = 50, halign = :left)

    
    # Playback loop state
    is_playing = Observable(false)
    
    on(play_btn.clicks) do _
        is_playing[] = !is_playing[]
        play_btn.label[] = is_playing[] ? "Pause" : "Play"
    end
    
    on(reset_btn.clicks) do _
        is_playing[] = false
        play_btn.label[] = "Play"
        set_close_to!(t_slider, sol.t[1])
    end
    
    # Asynchronous background loop for updating animation frames
    Base.@async begin
        fps = 30
        dt_real = 1.0 / fps
        has_opened = false
        
        while true
            # Monitor window open/close states safely
            if !has_opened
                if events(fig.scene).window_open[]
                    has_opened = true
                end
            else
                # Terminate loop when the window is closed to avoid memory/task leaks
                if !events(fig.scene).window_open[]
                    break
                end
            end
            
            sleep(dt_real)
            
            if is_playing[]
                t_curr = t_slider.value[]
                speed = speed_slider.value[]
                
                # Base speed scale: entire simulation plays in 10 seconds at 1x speed
                total_duration = sol.t[end] - sol.t[1]
                sim_dt = (total_duration / 10.0) * speed * dt_real
                
                t_next = t_curr + sim_dt
                if t_next >= sol.t[end]
                    set_close_to!(t_slider, sol.t[end])
                    is_playing[] = false
                    play_btn.label[] = "Play"
                else
                    set_close_to!(t_slider, t_next)
                end
            end
        end
    end
    
    return fig
end

"""
    save_orbit_animation(sol, bodies::Vector{Body}, filename::String; fps=30, duration=10.0, trail_duration=nothing, title="2D Gravitational Orbits", xlims=nothing, ylims=nothing)

Renders the orbit simulation and saves it to a GIF/MP4 file.

# Arguments
- `sol`: The ODESolution from `simulate_orbits`.
- `bodies`: The vector of `Body` instances.
- `filename`: Target output path (e.g. "orbits.gif" or "orbits.mp4").
- `fps`: Frame rate of the saved animation (default: 30).
- `duration`: Real-world playback duration of the video in seconds (default: 10.0).
- `trail_duration`: How long the orbital tail should be in simulation time units.
- `title`: Plot title.
- `xlims`, `ylims`: Fixed axis limits to override automatic scaling.
"""
function save_orbit_animation(
    sol,
    bodies::Vector{Body},
    filename::String;
    fps::Int=30,
    duration::Real=10.0,
    trail_duration::Union{Nothing, Real}=nothing,
    title::String="2D Gravitational Orbits",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing
)
    N = length(bodies)
    
    # Precompute default limits
    xs = Float64[]
    ys = Float64[]
    for t in range(sol.t[1], sol.t[end], length=300)
        state = sol(t)
        for i in 1:N
            push!(xs, state[4*(i-1) + 1])
            push!(ys, state[4*(i-1) + 2])
        end
    end
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)
    dx = xmax - xmin
    dy = ymax - ymin
    if dx ≈ 0.0; dx = 1.0; end
    if dy ≈ 0.0; dy = 1.0; end
    padding = 0.1
    default_xlims = (xmin - padding*dx, xmax + padding*dx)
    default_ylims = (ymin - padding*dy, ymax + padding*dy)

    # Setup Figure and Axis
    fig = Figure(size = (800, 650))
    ax = Axis(
        fig[1, 1],
        aspect = DataAspect(),
        title = title,
        xlabel = "x",
        ylabel = "y"
    )
    
    if !isnothing(xlims)
        xlims!(ax, xlims...)
    else
        xlims!(ax, default_xlims...)
    end
    
    if !isnothing(ylims)
        ylims!(ax, ylims...)
    else
        ylims!(ax, default_ylims...)
    end

    t_obs = Observable(sol.t[1])
    
    pos_observables = [
        lift(t -> get_position_at_t(sol, t, i), t_obs) for i in 1:N
    ]
    trail_observables = [
        lift(t -> get_trail_at_t(sol, t, i, trail_duration), t_obs) for i in 1:N
    ]
    
    for i in 1:N
        color_val = bodies[i].color
        trail_color = color_val isa Symbol ? (color_val, 0.6) : color_val
        
        lines!(ax, trail_observables[i], color = trail_color, linewidth = 2)
        scatter!(ax, pos_observables[i], color = color_val, markersize = bodies[i].size, marker = :circle, strokecolor = :white, strokewidth = 1.0)
        
        if !isempty(bodies[i].name)
            text!(ax, pos_observables[i], text = bodies[i].name, align = (:center, :bottom), offset = (0, 8), fontsize = 12, color = :black)
        end
    end
    
    # Label showing the current simulation time
    Label(fig[2, 1], @lift(string("Simulation Time: ", round($t_obs, digits=2))), halign = :center)
    
    # Render loop
    t_min, t_max = sol.t[1], sol.t[end]
    num_frames = Int(round(fps * duration))
    t_vals = range(t_min, t_max, length = num_frames)
    
    record(fig, filename, t_vals; fps = fps) do t
        t_obs[] = t
    end
    
    return nothing
end
