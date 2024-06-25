function sumSquareDifference(N::Integer)
    sumSquare = sum((1 : N).^2)
    squareSum = N^2 * (N + 1)^2 / 4

    return squareSum - sumSquare
end
