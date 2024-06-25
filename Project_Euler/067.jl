triangleString = read("InputFiles/067_triangle.txt", String)
triangle = Vector{Vector{Int}}()
for line ∈ split(triangleString, '\n')[1 : end - 1]
    vals = parse.(Int, split(line, ' '))
    push!(triangle, vals)
end

function maxPathSum(tri)
    n = length(tri)
    dp = deepcopy(tri)

    for row ∈ n - 1: -1 : 1
        for col ∈ 1 : row
            dp[row][col] += max(dp[row + 1][col], dp[row + 1][col + 1])
        end
    end

    return dp[1][1]
end

maxPathSum(triangle)
