% Generalized Gregory quadrature for computing high order hypergeometric pFq over a grid
%
% Author: Caleb Jacobs
% DLM: October 14, 2024

function [z, f, h] = pFq(a, b, varargin)
    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'r', 2.49);     % Width/height of grid
    addParameter(p, 'n', 60);       % Number of nodes for width/height
    addParameter(p, 'np', 3);       % Number of padding layers
    addParameter(p, 'Tr', 0.6);     % Taylor expansion radius
    addParameter(p, 'corrR', .5);   % z = 1 radius to correct
    addParameter(p, 'innerR', .6);  % Inner annulus radii for correction nodes
    addParameter(p, 'outerR', .8);  % Outer annulus radii for correction nodes
    addParameter(p, 'z1N', 70);     % Number of z = 1 expansion terms
    parse(p, varargin{:});

    r  = p.Results.r;
    n  = p.Results.n;
    np = p.Results.np;
    Tr = p.Results.Tr;
    corrR  = p.Results.corrR;
    innerR = p.Results.innerR;
    outerR = p.Results.outerR;
    z1N    = p.Results.z1N;

    o = min(length(a), length(b)); % Order of pFq

    g = Grid(n, r, Tr, np, o); % Initial grid

    aIdx = length(a); % Index for a
    bIdx = length(b); % Index for b

    branch = false; % Initialize extra branch cut correction
    if length(a) == length(b)
        f = exp(g.z);                  % 0F0 base function case
        h = f;
    elseif length(a) == length(b) + 1
        f = oneMinusZ_alpha(g.z, -a(end), false); % 1F0 base function case
        h = oneMinusZ_alpha(g.z, -a(end), true);  % 1F0 base function case alternate branch

        wa = [1];        % Set initial singular weight

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
                (gamfun(b(bIdx)) / ...
                 gamfun(a(aIdx)) / ...
                 gamfun(b(bIdx) - a(aIdx)));
        
        % Compute next set of pFq values
        if ~branch
            f = Gamma .* (D0 * f + D1 * h); % Compute values for the next layer
            h = f;
        else
            fTmp = Gamma .* (D0 * f + D1 * h); % Compute values for the next layer
            h = Gamma .* (D2 * f + D3 * h);    % Compute values for the next layer (alternate branch)
            f = fTmp;

            %% Use z = 1 expansion to correct
            zm1  = abs(g.z - 1);                            % Shift all z values by 1
            cIdx = zm1 <= corrR;                            % Corrected node indices

            zIdx = (innerR <= zm1) & (zm1 < outerR);        % Stencil nodes
            zCirc = g.z(zIdx);
            fCirc = f(zIdx);

            [wa, wb] = getZ1ExpansionWeights(a(aIdx : end), b(bIdx : end), wa, zCirc, fCirc, z1N);

            f(cIdx) = z1PFQ(a(aIdx : end), b(bIdx : end), wa, wb, g.z(cIdx), false);
            h(cIdx) = z1PFQ(a(aIdx : end), b(bIdx : end), wa, wb, g.z(cIdx), true);
        end

        f(g.c) = 1 + 0i; % Correct origin value

        aIdx = aIdx - 1; % Move to the next layer in a
        bIdx = bIdx - 1; % Move to the next layer in b
    end

    z = g.z;
end
