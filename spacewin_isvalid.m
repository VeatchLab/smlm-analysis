function valid = spacewin_isvalid(W)
types = {'polygon'};
valid = false;

if isstruct(W) && isfield(W,'type')
    valid = true;
    if ~any(validatestring(W.type, types))
        warning('spatial window type not known');
    end
end
