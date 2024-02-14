function  complexPlot3D(z, f)
% Create a complex abs-phase plot
Z = reshape(z,[],sqrt(length(z)));
F = reshape(f,[],sqrt(length(z)));

X = real(Z);
Y = imag(Z);
C = angle(F);
C(C < 0) = C(C < 0) + 2*pi;
F = abs(F);

surf(X, Y, F, C)
colormap('hsv')
clim([0, 2*pi])
end

