        %% getIdx function
        function zIdx = getIdx(z, g)
            if abs(real(z)) > g.r || abs(imag(z)) > g.r
                error('Argument must be inside of the grid');
            end

            Nx = round(real(z) / g.h);
            Ny = round(imag(z) / g.h);

            zIdx = g.c + Nx * g.dx + Ny * g.dy;

            if abs(z - g.z(zIdx)) ~= 0
                error('Argument must be on the complex grid');
            end
        end
