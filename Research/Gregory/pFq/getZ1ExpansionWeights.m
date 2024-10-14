function [Ak ,Bk] = getZ1ExpansionWeights(a, b, wa, z, f, n)
%GETZ1EXPANSIONWEIGHTS
%   Compute singular and regular expansion weights, Ak and Bk, given the
%   old singular weights and accurate functions values f at z. 
%   The a and b vectors are the pFq parameter vectors.

c = a(1);                       % Current a coefficient
d = b(1);                       % Current b coefficient
ga = sum(b) - sum(a) + c - d;   % Branch exponent

Ak = zeros(n,1);
for k = 0 : n - 1
    s = 0;
    for j = 0 : min(k, length(wa) - 1)
        s = s + (-1)^j * wa(j + 1) * binom(d + ga + k - 1, d + ga + j - 1) * ...
                                     binom(d - c + k - j - 1, - ga - j - 1);
    end

    Ak(k + 1) = pi * csc(pi * (c - ga - d)) * (d - c) * binom(d - 1, c - 1) * s;
end

sing = zeros(length(z),1);
for i = 1 : length(z)
    zm1 = oneMinusZ_alpha(z(i), ga - c + d, false);

    S = sum(Ak .* ((1 - z(i)).^transpose([0 : length(Ak) - 1])));

    sing(i) = zm1 * S;
end

A = (1 - z).^[0 : n - 1];
Bk = A \ (f - sing);
