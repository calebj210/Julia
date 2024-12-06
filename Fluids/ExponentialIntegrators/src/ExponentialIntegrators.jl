module ExponentialIntegrators

using DifferentialEquations,
      SciMLOperators,
      FFTW

using GLMakie

export kuramoto,
       plot_kuramoto,
       animate_kuramoto
include("kuramoto_sivashinsky.jl")

export practice, 
       lorenz
include("practice.jl")

export convolve
include("convolve.jl")
end
