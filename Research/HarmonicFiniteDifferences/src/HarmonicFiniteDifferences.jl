module HarmonicFiniteDifferences

include("centered_stencils.jl")
include("onesided_stencils.jl")

export sphericalharmonic,
       centered_nodes,
       onesided_nodes,
       onesided_sparse_structure,
       vand

end
