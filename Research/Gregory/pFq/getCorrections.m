% Compute all possible end/corner corrections over a grid `g`.
function result = getCorrections(g, alpha, beta)
    bpci = getBasepointCorrectionIndices(g); % Basepoint correction indices

    % Basepoint corrections
    bpr = getCorrectionWeights(g, bpci, alpha,  1 + 0i); % Basepoint right
    bpu = getCorrectionWeights(g, bpci, alpha,  0 + 1i); % Basepoint up
    bpl = getCorrectionWeights(g, bpci, alpha, -1 + 0i); % Basepoint left
    bpd = getCorrectionWeights(g, bpci, alpha,  0 - 1i); % Basepoint down

    % Corner corrections
    crr = getCorrectionWeights(g, bpci, 0,  1 + 0i); % Corner right
    cru = getCorrectionWeights(g, bpci, 0,  0 + 1i); % Corner up
    crl = getCorrectionWeights(g, bpci, 0, -1 + 0i); % Corner left
    crd = getCorrectionWeights(g, bpci, 0,  0 - 1i); % Corner down

    % Endpoint corrections
    epr = getCorrectionWeights(g, bpci, beta,  1 + 0i); % Endpoint right
    epu = getCorrectionWeights(g, bpci, beta,  0 + 1i); % Endpoint up
    epl = getCorrectionWeights(g, bpci, beta, -1 + 0i); % Endpoint left
    epd = getCorrectionWeights(g, bpci, beta,  0 - 1i); % Endpoint down

    result = Corrections(bpci, bpu, bpd, bpl, bpr, cru, crd, crl, crr, epu, epd, epl, epr);
end