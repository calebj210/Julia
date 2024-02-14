% Generate differentiation matrix for ∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du over a grid of radius `r`.
function [D0, D1, D2, D3] = getDiffMat(n, r, alpha, beta, ir, np, nl, branch)
g = Grid(n, r, ir, np, nl); % Generate grid

M = length(g.i) + length(g.e);
N = length(g.z);

if ~branch
    row0 = [];
    row1 = [];

    col0 = [];
    col1 = [];

    wgt0 = [];
    wgt1 = [];
else
    % Initialize row vectors
    row0 = [];
    row1 = [];
    row2 = [];
    row3 = [];

    % Initialize column vectors
    col0 = [];
    col1 = [];
    col2 = [];
    col3 = [];

    % Initialize weight vectors
    wgt0 = [];
    wgt1 = [];
    wgt2 = [];
    wgt3 = [];
end


% Populate internal weights using Taylor expansion approximation
% A = lu(getVand(g.ib, g)); % Compute vandermonde of internal boundary nodes
A = getVand(g.ib, g); % Compute vandermonde of internal boundary nodes

% Populate external weights using generalized Gregory quadrature
c = getCorrections(g, alpha, beta); % End/corner corrections

i = 1;

for j = 1:length(g.ie)
    t = g.ie(j);

    if t == 'p'
        continue;
    end

    if t == 'i'
        row0 = [row0, g.ib];
        col0 = [col0, repmat(i, 1, length(g.ib))];
        wgt0 = [wgt0, getInternalWeights(j, A, g, alpha, beta)];

        i = i + 1;
    elseif t == 'e'
        [idx0, wgts0, idx1, wgts1] = getExternalWeights(j, c, g, alpha, beta, false);

        row0 = [row0, idx0];
        col0 = [col0, repmat(i, 1, length(idx0))];
        wgt0 = [wgt0, wgts0];

        row1 = [row1, idx1];
        col1 = [col1, repmat(i, 1, length(idx1))];
        wgt1 = [wgt1, wgts1];

        if branch && real(g.z(i)) >= 1
            [idx2, wgts2, idx3, wgts3] = getExternalWeights(j, c, g, alpha, beta, true);

            row2 = [row2, idx2];
            col2 = [col2, repmat(i, 1, length(idx2))];
            wgt2 = [wgt2, wgts2];

            row3 = [row3, idx3];
            col3 = [col3, repmat(i, 1, length(idx3))];
            wgt3 = [wgt3, wgts3];
        end

        i = i + 1;
    end
end


% Construct sparse arrays
D0 = sparse(col0, row0, wgt0, M, N);
D1 = sparse(col1, row1, wgt1, M, N);
if branch
    D2 = sparse(col2, row2, wgt2, M, N);
    D3 = sparse(col3, row3, wgt3, M, N);
else
    D2 = [];
    D3 = [];
end

