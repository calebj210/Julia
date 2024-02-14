% Generalized, grid-based Gregory quadrature for computing hypergeometric pFq
%
% Author: Caleb Jacobs
% DLM: January 22, 2023

classdef Corrections
    properties
        bpci

        bpu
        bpd
        bpl
        bpr

        cru
        crd
        crl
        crr

        epu
        epd
        epl
        epr
    end

    methods
        %% getGrid function
        function corr = Corrections(bpci, bpu, bpd, bpl, bpr, cru, crd, crl, crr, epu, epd, epl, epr)
            corr.bpci=bpci;
            corr.bpu=bpu;
            corr.bpd=bpd;
            corr.bpl=bpl;
            corr.bpr=bpr;
            corr.cru=cru;
            corr.crd=crd;
            corr.crl=crl;
            corr.crr=crr;
            corr.epu=epu;
            corr.epd=epd;
            corr.epl=epl;
            corr.epr=epr;
        end
    end
end












