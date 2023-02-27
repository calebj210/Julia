# Find solutions to 3×3×3 block puzzle
#
# Author: Caleb Jacobs
# DLM: 27/02/2023

# Puzze data structure
mutable struct Chain
    const n::Int                # Number of links
    const link::Vector          # link data vector
    dirs::Array{Int64, 1}       # Direction at each link
    idx::Int                    # Current link
    pos::CartesianIndex{3}      # Current position
    state::Array{Int64, 3}      # Solution state tensor
end

# Check if chain is valid
function validChain(chain::Vector)::Bool
    n = length(chain)

    return (sum(chain) + n) == 27
end

# Toggle state under current position
function toggleState(C::Chain)::Bool
    for i ∈ 1 : 3
        if C.pos[i] < 1 || C.pos[i] > 3
            return false
        end
    end
    
    return C.state[C.pos] = C.state[C.pos] ⊻ true

end

# Get state under current positon
function getState(C::Chain)::Bool
    for i ∈ 1 : 3
        if C.pos[i] < 1 || C.pos[i] > 3
            return false
        end
    end

    return C.state[C.pos]
end

# Initialize chain
function initChain(chain)::Chain
    C = Chain(
            length(chain), 
            chain, 
            zeros(Int, length(chain)),
            1, 
            CartesianIndex(1,1,1), 
            BitArray(undef, 3, 3, 3))
    
    toggleState(C)
    C.state = C.state .⊻ true

    return C
end

# Try to walk in direction (true = success, false = failed)
function walk(C::Chain; rev = false)::Bool
    dir = C.dirs[C.idx]     # Get current walking direction
    
    if rev
        dir *= -1
    end

    # Translate direction into meaningful data
    if dir == 1
        v = CartesianIndex( 1, 0, 0)
    elseif dir == -1
        v = CartesianIndex(-1, 0, 0)
    elseif dir == 2
        v = CartesianIndex( 0, 1, 0)
    elseif dir == -2
        v = CartesianIndex( 0,-1, 0)
    elseif dir == 3
        v = CartesianIndex( 0, 0, 1)
    elseif dir == -3
        v = CartesianIndex( 0, 0,-1)
    end

    clear = true            # Is path currently clear

    # Check paths
    i = 0
    while clear && i < C.n
        C.pos += v
        toggleState(C)

        clear = getState(C)

        i += 1
    end

    # Check if path was walked successfully or return state
    if !clear
        while i >= 0
           toggleState(C)
           C.pos -= v
           i -= 1
        end
    end

    return clear
end

# Recursive chain solver
function solve(C::Chain)
    println(C.idx)

    C.idx += 1          # Move to next link

    # Check if we are at the end of the link
    if C.idx >= C.n
        return true
    end

    # Try moving in each direction
    for i ∈ 1 : 6
        C.dirs[C.idx] = i
        if walk(C)
            return solve(C)
        end
    end

    C.idx -= 1              # Couldn't move from here

    return false
end

# Start solving chain
function solve(chain::Array{Int64, 1})
    # Check if chain is valid
    if !validChain(chain)
        print("Invalid chain!")
        return 1
    end

    # Initialize solution chain
    C = initChain(chain)

    solved = false      # Assume there is no solution
    for i ∈ 1 : 4

        # Try to find solution
        if i == 2
            toggleState(C)
            C.pos = CartesianIndex(2,1,1)
            toggleState(C)
        elseif i == 3
            toggleState(C)
            C.pos = CartesianIndex(2,2,1)
            toggleState(C)
        elseif i == 4
            toggleState(C)
            C.pos = CartesianIndex(2,2,2)
            toggleState(C)
        end

        if solve(C)
            solved = true
            break
        end
    end
        
    return solved
end


function driver()
    chain = [1,0,0,1,0,1,0,0,1,1,1,0,0,1,1,1,1]

    solve(chain)
end
