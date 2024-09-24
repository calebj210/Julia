module ComplexVisuals

export complex_grid, 
       complex_phase_plot, 
       complex_surface_plot,
       complex_real_imag_plot

include("complex_builder.jl")
include("complex_2d_plots.jl")
include("complex_3d_plots.jl")

end
