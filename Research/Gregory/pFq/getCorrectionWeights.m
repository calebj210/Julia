function weights = getCorrectionWeights(g, idx, alpha, dir)
    N = length(idx);                                         % Number of correction nodes

    A = getVand(idx, g);                                     % Left hand side

    b = -transpose(arrayfun(@(k) zeta(-alpha - k) * (dir * g.h)^k, 0:N-1));  % Right hand side

    weights = A \ b;                                         % Solve for weights

    weights = weights * g.h^(1 + alpha);                     % Scale weights by grid spacing factor
end