function b = binom(n, k)
%BINOM
%   Compute the binomial coefficient nCk
%   for any complex arguments n and k.
b = gamfun(n + 1) / gamfun(k + 1) / gamfun(n - k + 1);
