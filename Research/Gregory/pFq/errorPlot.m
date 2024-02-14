function p = errorPlot(z, f, tru)
%Log Absolute error plot
Z = reshape(z, [], sqrt(length(z)));
F = reshape(f - tru, [], sqrt(length(f)));
C = log10(abs(F));
C(isinf(C)) = -16;

p = surf(real(Z), imag(Z), C, ...
         'EdgeColor', 'none');
colormap viridis
colorbar
clim([-16,1])
view(2)
pbaspect([1 1 1])
end