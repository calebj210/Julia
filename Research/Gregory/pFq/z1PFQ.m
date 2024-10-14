function f = z1PFQ(a, b, wa, wb, z, branch)
%Z1PFQ
%   Compute pFq(a;b;z) near z = 1 given expansion weights wa and wb.
al = sum(b) - sum(a);

f = zeros(length(z), 1);
for i = 1 : length(z)
    ga = oneMinusZ_alpha(z(i), al, branch);
    zi = (1 - z(i)).^transpose(0 : length(wa) - 1);

    f(i) = sum(ga * wa .* zi) + sum(wb .* zi);
end
