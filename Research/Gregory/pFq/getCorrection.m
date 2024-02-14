function result = getCorrection(c, dir, type)
switch type
    case 'cr'
        if isequal(dir, 0 + 1i)
            result = c.cru;
        elseif isequal(dir, 0 - 1i)
            result = c.crd;
        elseif isequal(dir, -1 + 0i)
            result = c.crl;
        elseif isequal(dir, 1 + 0i)
            result = c.crr;
        else
            error('Invalid direction of grid integration');
        end
    case 'bp'
        if isequal(dir, 0 + 1i)
            result = c.bpu;
        elseif isequal(dir, 0 - 1i)
            result = c.bpd;
        elseif isequal(dir, -1 + 0i)
            result = c.bpl;
        elseif isequal(dir, 1 + 0i)
            result = c.bpr;
        else
            error('Invalid direction of grid integration');
        end
    case 'ep'
        if isequal(dir, 0 + 1i)
            result = c.epu;
        elseif isequal(dir, 0 - 1i)
            result = c.epd;
        elseif isequal(dir, -1 + 0i)
            result = c.epl;
        elseif isequal(dir, 1 + 0i)
            result = c.epr;
        else
            error('Invalid direction of grid integration');
        end
    otherwise
        error('Invalid point type, needs a valid point type of bp, cr, or ep.');
end

result = transpose(result);
end