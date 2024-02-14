% Generalized Gregory quadrature for computing high order hypergeometric pFq over a grid
%
% Author: Caleb Jacobs
% DLM: December 2, 2023

function [z, f, h] = pFq(a, b, varargin)
    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'r', 1.99);
    addParameter(p, 'n', 40);
    addParameter(p, 'np', 3);
    addParameter(p, 'Tr', 0.6);
    parse(p, varargin{:});

    r  = p.Results.r;
    n  = p.Results.n;
    np = p.Results.np;
    Tr = p.Results.Tr;

    o = min(length(a), length(b)); % Order of pFq

    g = Grid(n, r, Tr, np, o); % Initial grid
    % g = getGrid(n, r, 'ir', Tr, 'np', np, 'nl', o); % Initial grid

    aIdx = length(a); % Index for a
    bIdx = length(b); % Index for b

    branch = false; % Initialize extra branch cut correction
    if length(a) == length(b)
        f = exp(g.z);                  % 0F0 base function case
        h = f;
    elseif length(a) == length(b) + 1
        f = oneMinusZ_alpha(g.z, -a(end), false); % 1F0 base function case
        h = oneMinusZ_alpha(g.z, -a(end), true);  % 1F0 base function case alternate branch

        branch = true;

        aIdx = aIdx - 1; % Move a index to the next layer
    elseif length(a) == length(b) - 1
        f = hypergeom([],b(bIdx), g.z); % 0F1 base function case
        h = f;

        bIdx = bIdx - 1; % Move b index to the next layer
    else
        error('Non-supported pFq order');
    end

    for nl = o : -1 : 1
        alpha = a(aIdx) - 1; % Base point singularity order
        beta = b(bIdx) - a(aIdx) - 1; % Evaluation point singularity order

        [D0, D1, D2, D3] = getDiffMat(n, r, alpha, beta, Tr, np, nl, branch);

        g = Grid(n, r, Tr, np, nl - 1); % Get next computation grid

        Gamma = g.z.^(1 - b(bIdx)) * ... % Front coefficient
                (gamma(b(bIdx)) / ...
                 gamma(a(aIdx)) / ...
                 gamma(b(bIdx) - a(aIdx)));

        aIdx = aIdx - 1; % Move to the next layer in a
        bIdx = bIdx - 1; % Move to the next layer in b
        
        if ~branch
            f = Gamma .* (D0 * f + D1 * h); % Compute values for the next layer
            h = f;
        else
            fTmp = Gamma .* (D0 * f + D1 * h); % Compute values for the next layer
            h = Gamma .* (D2 * f + D3 * h); % Compute values for the next layer
            f = fTmp;
        end

        f(g.c) = 1 + 0i;
    end

    z = g.z;
end
