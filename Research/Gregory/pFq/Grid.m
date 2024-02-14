% Constructors and routines for working with complex grids for quadrature.
%
% Author: Caleb Jacobs
% DLM: January 22, 2023

%% Complex grid for use with grid-based quadratures
classdef Grid
    properties
        z % Grid points
        i % Internal indices
        e % External indices
        ie % Internal-external-padding vector
        ib % Internal node boundary
        dx % Index spacing in x
        dy % Index spacing in y
        c % Origin index
        h % Grid spacing
        r % Radius of square grid
        np % Number of padded nodes
        nl % Number of padding layers
        T % Index radius of Taylor expansion
    end

    methods
        %% getGrid function
        function grid = Grid(n, r, ir, np, nl)
            x = linspace(0, r, n + 1);
            y = linspace(0, r, n + 1);
            h = abs(x(2) - x(1));

            pad = r + h * (1:nl*np);
            xp = [-pad(end:-1:1), -x(end:-1:2), x, pad];
            yp = [-pad(end:-1:1), -y(end:-1:2), y, pad];

            [X, Y] = meshgrid(xp, yp);
            gr = X + 1i * Y;

            %dx = stride(grid, 2);
            %dy = stride(grid, 1);
            dx = size(gr,1);
            dy = 1;

            c = 1 + (np * nl + n) * (dx + dy);
            T = round(ir / h);

            if nl==0
                pr = r;
            else
                pr = r + h * (nl - 1) * np;
            end


            z = gr(:);
            i = [];
            e = [];
            ie = [];
            ib = [];
            p =  [];

            for idx = 1:length(z)
                if abs(z(idx)) <= ir
                    i = [i, idx];
                    ie = [ie, 'i'];

                    if abs(z(idx)) > ir - h
                        ib = [ib, idx];
                    end
                elseif abs(real(z(idx))) <= pr && abs(imag(z(idx))) <= pr
                    e = [e, idx];
                    ie = [ie, 'e'];
                else
                    ie = [ie, 'p'];
                end
            end

            [~, order] = sort(angle(z(ib)));
            ib = ib(order);

            grid.z = z;
            grid.i = i;
            grid.e = e;
            grid.ie = ie;
            grid.ib = ib;
            grid.dx = dx;
            grid.dy = dy;
            grid.c = c;
            grid.h = h;
            grid.r = r;
            grid.np = np;
            grid.nl = nl;
            grid.T = T;

            %grid = Grid(z, i, e, ie, ib, dx, dy, c, h, r, np, nl, T);
        end


%         %% getGrid function
%         function grid = getGrid(n, r, ir, np, nl)
%             x = linspace(0, r, n + 1);
%             y = linspace(0, r, n + 1);
%             h = abs(x(2) - x(1));
% 
%             pad = r + h * (1:nl*np);
%             xp = [-pad(end:-1:1), -x(end:-1:2), x, pad];
%             yp = [-pad(end:-1:1), -y(end:-1:2), y, pad];
% 
%             [X, Y] = meshgrid(xp, yp);
%             grid = X + 1i * Y;
% 
%             dx = stride(grid, 2);
%             dy = stride(grid, 1);
%             c = 1 + (np * nl + n) * (dx + dy);
%             T = round(ir / h);
% 
%             %pr = nl == 0 ? r : r + h * (nl - 1) * np;
% 
%             if nl==0
%                 pr = r;
%             else
%                 pr = r + h * (nl - 1) * np;
%             end
% 
% 
%             z = grid(:);
%             i = [];
%             e = [];
%             ie = [];
%             ib = [];
% 
%             for idx = 1:length(z)
%                 if abs(z(idx)) <= ir
%                     i = [i, idx];
%                     ie = [ie, 'i'];
% 
%                     if abs(z(idx)) > ir - h
%                         ib = [ib, idx];
%                     end
%                 elseif abs(real(z(idx))) <= pr && abs(imag(z(idx))) <= pr
%                     e = [e, idx];
%                     ie = [ie, 'e'];
%                 else
%                     ie = [ie, 'p'];
%                 end
%             end
% 
%             [~, order] = sort(angle(z(ib)));
%             ib = ib(order);
% 
%             grid = Grid(z, i, e, ie, ib, dx, dy, c, h, r, np, nl, T);
%         end

    end

end
