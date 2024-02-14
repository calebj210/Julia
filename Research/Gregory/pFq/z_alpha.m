% Compute roots given a power α and a branch cut rotation of θ.
function result = z_alpha(Z, alpha, theta)
result = zeros(size(Z));
for i = 1 : length(Z)
    z = Z(i);
    
    if theta >= 0
        if angle(z) >= theta - pi
            result(i) = z^alpha;
        else
            result(i) = z^alpha .* cispi(2 * mod(alpha, 1));
        end
    else
        if angle(z) < theta + pi
            result(i) = z^alpha;
        else
            result(i) = z^alpha .* cispi(-2 * mod(alpha, 1));
        end
    end
end
