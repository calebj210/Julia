module ComplexVisuals

export complex_color_wheel,
       complex_color_wheel!,
       complex_theme

export ComplexGrid,
       complex_square_grid

export Phase,
       phase,
       phase!

export ComplexSurface,
       complexsurface,
       complexsurface!

using MakieCore, Makie

include("complex_grids.jl")
include("complex_range.jl")
include("color_wheel.jl")
include("complex_plots_2d.jl")
include("complex_plots_3d.jl")
end
