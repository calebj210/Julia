% Get extra branch rotation correction
function result = theta_gamma(z, ze, gamma)
if imag(z) == 0
    if imag(ze) >= 0
        result = cispi(2 * mod(gamma, 1));
    else
        result = 1;
    end
elseif sgn(imag(z)) == sgn(imag(ze))
    result = 1;
else
    result = cispi(-2 * sgn(imag(z)) * mod(gamma, 1));
end
end