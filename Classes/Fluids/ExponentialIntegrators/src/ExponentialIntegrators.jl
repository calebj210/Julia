module ExponentialIntegrators

using DifferentialEquations,
      SciMLOperators,
      FFTW

using CairoMakie
using LinearAlgebra

export kuramoto,
       plot_kuramoto,
       animate_kuramoto
include("kuramoto_sivashinsky.jl")

export burgers,
       burgers_heatmap,
       burgers_animation
include("burgers.jl")

export schrodinger,
       schrodinger_heatmap,
       schrodinger_animation
include("schrodinger.jl")

export KdV,
       plot_KdV,
       animate_KdV
include("KdV.jl")

export burgers_plot,
       schrodinger_plot,
       kuramoto_plot,
       KdV_plot
include("testing.jl")

export convolve
include("convolve.jl")
end
