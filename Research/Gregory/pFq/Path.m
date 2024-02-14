%% Path structure
classdef Path
    properties
        i
        p
        f
    end

    methods
        %% getPath function
        function path = Path(i, p, f)
            path.i = i;
            path.p = p;
            path.f = f;
        end

    end

end
