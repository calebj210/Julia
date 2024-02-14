% Compute appropriate branch cut rotation angles based on the path to `z`.
function [theta_alpha, theta_beta] = getBranchAngle(z, g, branch)
if ~branch
    if abs(real(z)) < abs(imag(z)) || (real(z) >= 0 && abs(imag(z)) >= 4 * g.np * g.h)
        % Right and left moving cuts
        if real(z) >= 0
            theta_alpha = sgn(imag(z)) * pi;
            theta_beta = 0;
        else
            theta_alpha = 0;
            theta_beta = sgn(imag(z)) * pi;
        end
    else
        % Up and down moving cuts
        if real(z) >= 0
            if imag(z) == 0
                theta_alpha = 0;
                theta_beta  = pi / 2;
            else
                theta_alpha = 0;
                theta_beta  = -sgn(imag(z)) * pi / 2;
            end
        else
            theta_alpha = sgn(imag(z)) * pi;
            theta_beta  = sgn(imag(z)) * pi / 2;
        end
    end
else
    theta_alpha = 0;
    theta_beta  = nsgn(imag(z));
end
end