function ind = spacewin_isinside(x,y,W)
% IND = SPACEWIN_ISINSIZE(X,Y,W)
%       IND(i) is true if the point X(i), Y(i) is in the spatial window W,
%       as specified by spacewin_isvalid(), and false otherwise.

if ~spacewin_isvalid(W)
    error('spacewin_isinside: invalid spatial window provided');
end

switch W.type
    case 'polygon'
        ind = inpolygon(x,y,W.x, W.y);
    case 'polyshape'
        [interior, onedge] = isinterior(W.p, x,y);
        ind = interior | onedge;
    case 'image'
        [r,c] = W.ref.worldToSubscript(x,y);
        i1 = ~isnan(r) & ~isnan(c);
        j = sub2ind(size(W.im), r(i1),c(i1));
        i2 = logical(W.im(j));
        ind = i1;
        ind(i1) = i2;
    otherwise
        error('spacewin_inside: spatial window of this type is not yet supported')
end
