function A = spacewin_area(W)
% A = SPACEWIN_AREA(W)      W should be a valid spatial window according to
%                           SPACEWIN_ISVALID(). A is the area of the
%                           spatial window, in appropriate units.

if ~spacewin_isvalid(W)
    error('spacewin_area: invalid spatial window')
end

switch W.type
    case 'polygon'
        A = polyarea(W.x, W.y);
    case 'polyshape'
        A = W.p.area();
    case 'image'
        A = sum(W.im(:))*W.ref.PixelExtentInWorldX*W.ref.PixelExtentInWorldY;
    otherwise
        error('spacewin_area: not implemented for this type of spatial window');
end