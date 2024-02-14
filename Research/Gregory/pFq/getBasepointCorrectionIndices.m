%% getBasepointCorrectionIndices function
function bpci = getBasepointCorrectionIndices(g)
    bpci = [];

    for i = g.c - g.dx * (g.np + 1) : g.c + g.dx * (g.np + 1)
        if (g.np - 1.5) * g.h <= abs(g.z(i)) && abs(g.z(i)) <= g.np * g.h
            bpci = [bpci, i];
        end
    end
    
    [~, I] = sort(angle(g.z(bpci)));
    bpci = bpci(I);
end