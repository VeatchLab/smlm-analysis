function spacewin = polyshape_to_spacewin(p, type, ref)

% TODO: check that the right kinds of variable are given
spacewin.type = type;
switch type
    case 'polyshape'
        spacewin.p = p;
        
    case 'image'
        spacewin.ref = ref;
        sz = ref.ImageSize;
        % correct, but very slow
%         [Xint,Yint] = meshgrid(1:sz(2), 1:sz(1));
%         [Xworld, Yworld] = ref.intrinsicToWorld(Xint,Yint);
%         BW = isinterior(p,Xworld(:),Yworld(:));
%         spacewin.im = reshape(BW,sz);
        BW = zeros(sz);
        holes = ishole(p);
        for i = 1:numel(holes)
            [xw,yw] = boundary(p,i);
            [xi,yi] = ref.worldToIntrinsic(xw,yw);
            mask = poly2mask(xi,yi,sz(1),sz(2));
            if holes(i)
                BW = BW & ~mask;
            else
                BW = BW | mask;
            end
        end
        spacewin.im = BW;
        
    case 'polygon'
        [xi,yi] = p.boundary(1); % use the first boundary
        if p.NumHoles ~= 0 || p.NumRegions ~= 1
            warning('polyshape cannot be converted to polygon, using first boundary (usually exterior)');
        end
        spacewin.x = xi;
        spacewin.y = yi;
    otherwise
        error('polyshape_to_spacewin: spatial window type not supported');
end