# Find solutions to 3×3×3 block puzzle
#
# Author: Caleb Jacobs
# DLM: 25/02/2023

mutable struct Chain
    const n::Int                # Number of links
    const link::Vector          # link data vector
    start::Array{Int64, 1}      # Start of chain
    dirs::Array{Int64, 1}       # Direction at each link
    idx::Int                    # Current link
    pos::Array{Int64, 1}        # Current position
    state::Array{Int64, 3}      # Solution state tensor
end

function validChain(chain::Vector)
    n = length(chain)

    return (sum(chain) + n) == 27
end

function initChain(chain)
    return Chain(
            length(chain), 
            chain, 
            [1,1,1]
            zeros(Int, length(chain))
            1, 
            zeros(Int, 3), 
            zeros(Int, 3, 3, 3))
end

function solve(chain)
    # Check if chain is valid
    if !validChain(chain)
        print("Invalid chain!")
        return 1
    end

    # Initialize solution chain
    C = initChain(chain)

    for i = 1 : 4
end


function driver()
    chain = [1,0,0,1,0,1,0,0,1,1,1,0,0,1,1,1,1]


    solve(Chain)
end
