function pth = getPath(zIdx, g, r, varargin)

p = inputParser;
addOptional(p, 'branch', false);
parse(p, varargin{:});
branch = p.Results.branch;

Nx = round(real(g.z(zIdx)) / g.h);
Ny = round(imag(g.z(zIdx)) / g.h);

if ~branch
    % Standard left-right U-contour
    if abs(Nx) < abs(Ny) || (Nx >= 0 && abs(Ny) >= 4 * g.np)
        h1 = g.c - g.dx * sgn(Nx) * (1:r - 1);
        c1 = h1(end) - g.dx * sgn(Nx);
        v = c1 + g.dy * sgn(Ny) * (1:abs(Ny) - 1);
        c2 = v(end) + g.dy * sgn(Ny);
        h2 = c2 + g.dx * sgn(Nx) * (1:r + abs(Nx) - 1);

        p1 = Path(g.c, h1, c1);
        p2 = Path(c1, v, c2);
        p3 = Path(c2, h2, zIdx);

        pth = [p1, p2, p3];

    % Shortened up-down U-contour for approaching branch cut
    elseif Nx >= 0
        if Ny == 0
            dir = -1;
        else
            dir = sgn(Ny);
        end

        v1 = g.c + g.dy * dir * (1:(4 * r - 1));

        c1 = v1(end) + g.dy * dir;

        h = c1 + g.dx * (1:abs(Nx) - 1);

        c2 = h(end) + g.dx;
        
        v2 = c2 - g.dy * dir * (1:(4 * r - abs(Ny) - 1));

        p1 = Path(g.c, v1, c1);
        p2 = Path(c1, h, c2);
        p3 = Path(c2, v2, zIdx);

        pth = [p1, p2, p3];
        
    % Standard up-down U-contour
    else
        v1 = g.c - g.dy * sgn(Ny) * (1:(3 * r - 1));
        c1 = v1(end) - g.dy * sgn(Ny);
        h = c1 + g.dx * sgn(Nx) * (1:(abs(Nx) - 1));
        c2 = h(end) + g.dx * sgn(Nx);
        v2 = c2 + g.dy * sgn(Ny) * (1:(3 * r + abs(Ny) - 1));

        p1 = Path(g.c, v1, c1);
        p2 = Path(c1, h, c2);
        p3 = Path(c2, v2, zIdx);

        pth = [p1, p2, p3];
    end
% Up-down U-contour that passes the branch cut
else
    v1 = g.c - g.dy * nsgn(Ny) * (1:(4 * r - 1));

    c1 = v1(end) - g.dy * nsgn(Ny);

    h = c1 + g.dx * sgn(Nx) * (1:(abs(Nx) - 1));
    
    c2 = h(end) + g.dx * sgn(Nx);

    v2 = c2 + g.dy * nsgn(Ny) * (1:(4 * r + abs(Ny) - 1));

    p1 = Path(g.c, v1, c1);
    p2 = Path(c1, h, c2);
    p3 = Path(c2, v2, zIdx);

    pth = [p1, p2, p3];

end
end
