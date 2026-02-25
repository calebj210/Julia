module HarmonicFiniteDifferences

include("stencil_solvers.jl")
include("centered_stencils.jl")
include("onesided_stencils.jl")
include("cardinalfunctions.jl")
include("visuals.jl")

export sphericalharmonic,
       circularharmonic,
       harmonic,
       centered_nodes,
       onesided_nodes,
       sparse_structure,
       vand,
       FD_weights,
       cardinalweights,
       cardinalfunction,
       onesided_2d,
       plot_weights,
       barplot,
       barplot!

end
