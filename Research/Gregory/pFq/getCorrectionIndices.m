%% getCorrectionIndices function
function correctionIndices = getCorrectionIndices(zIdx, bpci, g)
Nx = round(real(g.z(zIdx)) / g.h);
Ny = round(imag(g.z(zIdx)) / g.h);

correctionIndices = bpci + Nx * g.dx + Ny * g.dy;
end