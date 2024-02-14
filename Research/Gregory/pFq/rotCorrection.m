%% rotCorrection function
function rotatedIndices = rotCorrection(indices, dir)
r = length(indices) / 8;
rotatedIndices = circshift(indices, -2 * r * mod(dir, 4));
end