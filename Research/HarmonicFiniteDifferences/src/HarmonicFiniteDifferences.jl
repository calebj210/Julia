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
       cardinalweights,
       cardinalfunction,
       onesided_2d

end
