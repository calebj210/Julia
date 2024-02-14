function [pIdx, bIdx] = sortPath(idx, g, z, branch)
pIdx = [];
bIdx = [];

% No need for branching
if real(z) < 1
    pIdx = 1 : length(idx);
    return
end

% Decide whether the index is on one sheet or another
for i = 1:length(idx)
    v = idx(i);
    if ~branch
        if nsgn(imag(g.z(v))) == nsgn(imag(z))
            pIdx = [pIdx, i];
        else
            bIdx = [bIdx, i];
        end
    else
        if imag(z) == 0
            if nsgn(imag(g.z(v))) == 1
                pIdx = [pIdx, i];
            else
                bIdx = [bIdx, i];
            end
        elseif imag(z) > 0
            if nsgn(imag(g.z(v))) == -1
                pIdx = [pIdx, i];
            else
                bIdx = [bIdx, i];
            end
        else
            if nsgn(imag(g.z(v))) == 1
                pIdx = [pIdx, i];
            else
                bIdx = [bIdx, i];
            end
        end
    end
end
end



