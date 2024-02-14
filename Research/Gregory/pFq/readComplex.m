function f = readComplex(file)
%Read and parse a file of a complex numbers
A = readtable(file);
f = A.Var1 + A.Var2 * 1i;
end

