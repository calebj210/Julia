function f = z1Correction(a, b, wa, wb, z)
%Z1CORRECTION
%   Compute corrected values of f at z given singular coefficients wa and
%   regular coefficients wb.
N  = length(wa);                            % Number of singular coefficients
al = sum(b) - sum(a);                       % Alpha power

f = [];
for i = 1 : length(z)
    ga = oneMinusZ_alpha(z(i), al, false); % (1 - z)^alpha with correct branching
    tmp = sum((ga * wa + wb) .* (1 - z(i)).^transpose(0 : N - 1));
    f = [f; tmp];
end