module MineralProspecting

include("continuation.jl")
export 
    apply_operator,
    gradient_continue,
    two_layer_continue

include("test.jl")
export
    single_step_grad,
    single_step_layr,
    single_step_convergence_plot,
    multi_step_convergence_plot

end

