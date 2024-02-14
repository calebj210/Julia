        %% getReducedGridMap function
        function [iMap, eMap] = getReducedGridMap(g)
            iMap = zeros(size(g.i));
            eMap = zeros(size(g.e));

            iIdx = 1;
            eIdx = 1;
            tIdx = 1;

            while iIdx <= length(g.i) && eIdx < length(g.e)
                if g.i(iIdx) < g.e(eIdx)
                    iMap(iIdx) = tIdx;
                    iIdx = iIdx + 1;
                else
                    eMap(eIdx) = tIdx;
                    eIdx = eIdx + 1;
                end

                tIdx = tIdx + 1;
            end

            if iIdx <= length(g.i)
                iMap(iIdx:end) = tIdx:(length(g.i) + length(g.e));
            else
                eMap(eIdx:end) = tIdx:(length(g.i) + length(g.e));
            end
        end
