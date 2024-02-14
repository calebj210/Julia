function [rowfIdx, rowf, rowhIdx, rowh] = getExternalWeights(zIdx, c, g, alpha, beta, branch)

pth = getPath(zIdx, g, 2 * g.np, branch);
N = length(pth);
z = g.z(zIdx);

rowf = zeros(1, length(g.z));
rowh = zeros(1, length(g.z));

for n = 1:N
    p = pth(n);

    dir = sign(g.z(p.f) - g.z(p.i));

    iIdx = getCorrectionIndices(p.i, c.bpci, g);
    fIdx = getCorrectionIndices(p.f, c.bpci, g);

    % Fix indices if needed
    if n == N && real(z) > 1
        [ppIdx, bpIdx] = sortPath(p.p,  g, z, branch);
        [piIdx, biIdx] = sortPath(iIdx, g, z, branch);
        [pfIdx, bfIdx] = sortPath(fIdx, g, z, branch);
    else
        ppIdx = 1:length(p.p);
        bpIdx = [];
        piIdx = 1:length(iIdx);
        biIdx = [];
        pfIdx = 1:length(fIdx);
        bfIdx = [];
    end

    % Choose appropriate branch cuts
    [theta_alpha, theta_beta] = getBranchAngle(z, g, branch);

    % Trapezoidal weights
    alpha_t = transpose(z_alpha(g.z(p.p), alpha, theta_alpha));
    beta_t = transpose(z_alpha(z - g.z(p.p), beta, theta_beta));

    rowf(p.p(ppIdx)) = rowf(p.p(ppIdx)) + dir * g.h * alpha_t(ppIdx) .* beta_t(ppIdx);
    rowh(p.p(bpIdx)) = rowh(p.p(bpIdx)) + dir * g.h * alpha_t(bpIdx) .* beta_t(bpIdx);

    % Initial endpoint correction
    h = (n == 1) * z_alpha(dir, 1 + alpha, theta_alpha) + (n ~= 1) * dir;
    alpha_t = (n == 1) * ones(1, length(iIdx)) + (n ~= 1) * z_alpha(g.z(iIdx), alpha, theta_alpha);
    beta_t = transpose(z_alpha(z - g.z(iIdx), beta, theta_beta));

    if n == 1
        cond = 'bp';
    else
        cond = 'cr';
    end
    tmp = getCorrection(c, dir, cond);

    rowf(iIdx(piIdx)) = rowf(iIdx(piIdx)) + h * alpha_t(piIdx) .* beta_t(piIdx) .* tmp(piIdx);
    rowh(iIdx(biIdx)) = rowh(iIdx(biIdx)) + h * alpha_t(biIdx) .* beta_t(biIdx) .* tmp(biIdx);

    % Final endpoint correction
    h = (n == N) * z_alpha(dir, 1 + beta, theta_beta) + (n ~= N) * dir;
    alpha_t = transpose(z_alpha(g.z(fIdx), alpha, theta_alpha));
    beta_t = (n == N) * ones(1, length(fIdx)) + (n ~= N) * z_alpha(z - g.z(fIdx), beta, theta_beta);

    if n == N
        cond = 'ep';
    else
        cond = 'cr';
    end
    tmp = getCorrection(c, -dir, cond);
    
    rowf(fIdx(pfIdx)) = rowf(fIdx(pfIdx)) + h * alpha_t(pfIdx) .* beta_t(pfIdx) .* tmp(pfIdx);
    rowh(fIdx(bfIdx)) = rowh(fIdx(bfIdx)) + h * alpha_t(bfIdx) .* beta_t(bfIdx) .* tmp(bfIdx);
end

rowfIdx = find(rowf ~= 0);
rowf = rowf(rowfIdx);

rowhIdx = find(rowh ~= 0);
rowh = rowh(rowhIdx);
end