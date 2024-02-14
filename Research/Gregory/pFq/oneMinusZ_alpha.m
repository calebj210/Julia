function result = oneMinusZ_alpha(Z, alpha, branch)
% Compute roots given a power alpha and a branch cut rotation.
result = zeros(size(Z));
for i = 1 : length(Z)
    z = Z(i);
    if ~branch
        % if imag(z) == 0 && real(z) > 1
        %     result(i) = (1 - z)^alpha * cispi(2 * mod(alpha, 1));
        % else
            result(i) = (1 - z)^alpha;
        % end
    else
        if imag(z) == 0 && real(z) > 1
            result(i) = (1 - z)^alpha * cispi(-2 * mod(alpha, 1));
        elseif imag(z) > 0
            result(i) = (1 - z)^alpha * cispi(2 * mod(alpha, 1));
        else
            result(i) = (1 - z)^alpha * cispi(-2 * mod(alpha, 1));
        end
    end
end

