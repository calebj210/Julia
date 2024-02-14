% Compute internal weights for ∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du where `idx` are the indices of the internal nodes in the grid `g`.
function result = getInternalWeights(zIdx, A, g, alpha, beta)
z = g.z(zIdx); % Get z value at index
N = size(A, 1);

% Right hand side from Taylor expansion of integral
b = z.^(1 + alpha + beta + (0:N-1)) .* gamma(1 + alpha + (0:N-1)) .* gamma(1 + beta) ./ gamma(2 + alpha + beta + (0:N-1));

omega = transpose(A \ transpose(b));

result = omega;
end