% Generate vandermond type matrix from nodes in the grid `g` with indices `idx`.
function A = getVand(idx, g)
    A = zeros(length(idx), length(g.z(idx)));

    for i = 1:length(idx)
        A(:, i) = g.z(idx(i)).^(0:(length(idx)-1));
    end
end