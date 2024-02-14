function p = complexPlot(z, f)
%Create a complex phase portrait
Z = reshape(z, [], sqrt(length(z)));
F = reshape(f, [], sqrt(length(f)));
C = angle(F);
C(C < 0) = C(C < 0) + 2*pi;

p = surf(real(Z), imag(Z), C, ...
         'EdgeColor', 'none');
colormap hsv;
clim([0,2*pi])
view(2)
pbaspect([1 1 1])
end