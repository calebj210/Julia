function A = SPH(D,az,el,r)
% Input parameters
%    D  Highest SPH degree to include
%       Three parameters below: column vectors with evaluation points,
%       as obtained from x,y,z vectors by [az,el,r] = cart2sph(x,y,z)
%   az  azimuth:    polar angle in x,y-plane
%   el  elevation:  angle up/down from x,y-plane (not angle to z-axis)
%    r  radius:     distance of evaluation point from origin
% Output parameter
%   A   Array with one row for each evaluation point
%       (D+1)^2 columns with values for the successive SPH functions
%
% Starting with x,y,z-values, SPH returns in its leading columns
% (with r denoting sqrt(3)):
% d = 0:        1 ,
% d = 1:        y ,     z  ,           x       ,
% d = 2:     r x y,   r y z,    z^2-(x^2+y^2)/2,    r x z,  r(x^2-y^2)/2.

A = zeros(length(el(:)), (D+1)^2);  % Array to store the computed SPH values
for d = 0:D
    P = legendre(d,sin(el),'sch');
    A(:,d*(d+1)+1:(d+1)^2) = P';
end
for d = 0:D
    A(:,d^2+1:d*(d+1)) = fliplr(A(:,d*(d+1)+2:(d+1)^2).*sin((1:d).*az));
    A(:,d*(d+1)+1:(d+1)^2) = A(:,d*(d+1)+1:(d+1)^2) .*cos((0:d).*az);
    A(:,d^2+1:(d+1)^2) = A(:,d^2+1:(d+1)^2).*(r.^d);
end
