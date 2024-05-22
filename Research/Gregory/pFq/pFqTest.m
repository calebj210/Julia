function [z, f, h, tru] = pFqTest(a,b,r,n,np,Tr,sr,cr)
%Run a pFq test importing true solution
tru = readComplex('../Data/pfq.csv');
[z, f, h] = pFq(a, b, 'r', r, 'n', n, 'np', np, 'Tr', Tr, 'sr', sr, 'cr', cr);
errorPlot(z,f,tru);
end

