module MineralProspecting

include("continuation.jl")
export 
    apply_operator,
    gradient_continue,
    two_layer_continue

include("test.jl")
export
    continuation_compare,
    single_step_grad,
    single_step_two,
    single_step_grad_convergence_plot,
    single_step_two_convergence_plot

end

