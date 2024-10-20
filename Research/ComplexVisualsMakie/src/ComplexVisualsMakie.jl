module ComplexVisualsMakie

export CairoMakie,
       ComplexGrid,
       complex_grid,
       complex_phase_plot,
       complex_phase_plot!,
       complex_color_wheel,
       complex_color_wheel!,
       complex_surface_plot,
       complex_real_plot,
       complex_imag_plot

using CairoMakie
CairoMakie.activate!()

include("complex_grids.jl")
include("color_wheel.jl")
include("complex_plots_2d.jl")
include("complex_plots_3d.jl")
end
