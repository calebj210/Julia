function [wa ,wb] = getZ1ExpansionWeights(a, b, wa, z, f)
%GETZ1EXPANSIONWEIGHTS
%   Compute singular and regular expansion weights, wa and wb, given the
%   old singular weights and accurate functions values f at z. 
%   The a and b vectors are the pFq parameter vectors.
Nz = length(z);                 % Number of z values
Nw = length(wa);                % Number of old singular weights

c = a(1);                       % c parameter in expansion
d = b(1);                       % d parameter in expansion
al = sum(b) - sum(a) + c - d;   % Singular power of previous expansion

%% Compute phi (i.e. branch correction)
phi = 0;
for k = 0 : Nw - 1
    for l = 0 : 35
        phi = phi + ...
            wa(k + 1) * (-1)^l * (1 - z).^(k + l) * gamma(1 + k + l + al) / ...
            (factorial(l) * gamma(c - l) * gamma(1 + d - c + k + l + al));
    end
end
phi = 2i * (-1)^(d - c + al) * (1 - z).^(d - c + al) * ...
      sin(pi * al) * gamma(d) ./ z.^(d - 1) .* phi;


%% Compute weights
al = sum(b) - sum(a);   % Singular power of current expansion
A = (1 - z).^(0 : Nz - 1);                          % LHS of system
ba = phi ./ (1 - z).^al / (exp(2i * pi * al) - 1)   % Singular RHS
bb = f - phi / (exp(2i * pi * al) - 1)              % Regular RHS

wa = A \ ba;                                        % Singular weights
wb = A \ bb;                                        % Regular weights
end

