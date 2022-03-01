function correction = spatial_edge_correction(W, r)
% C = SPATIAL_EDGE_CORRECTION(W,R)      Isotropic spatial edge correction
%             for pair correlation functions or Ripley's K functions, for a
%             spatial window W (as specified for spacewin_isvalid()), at
%             distances R
%             Normalized to Area, so that it is order 1.

if strcmp(W.type, 'image')
    % Strategy is a bit different for images, for computational reasons
    wij = conv2_fft_norm_and_interp(W.im, W.ref);
else
    area = spacewin_area(W);
    wij = @(dx,dy) wij_poly(W,dx,dy)/area;
end

ntheta = 60;
thetas = pi*(1:ntheta)/ntheta; % only need 0 to pi because w_ij(theta + pi) = w_ij(theta)
correction = zeros(size(r));
for i = 1:numel(r)
    for j = 1:ntheta
        dx = cos(thetas(j))*r(i);
        dy = sin(thetas(j))*r(i);
        correction(i) = correction(i) + wij(dx, dy)/ntheta;
    end
end
end

function w = wij_poly(W, dx, dy)
switch W.type
    case 'polygon'
        mx = W.x;
        my = W.y;
        [x, y] = polybool('intersection', mx, my, mx + dx, my + dy);
        w = polyarea(x,y);
    case 'polyshape'
        p = W.p;
        pt = p.translate(dx, dy);
        pinter = intersect(p,pt);
        w = area(pinter);
end
end

function W = conv2_fft_norm_and_interp(I, ref)

sza = size(I, 1)*2 + 1;
szb = size(I, 2)*2 + 1;

C2 = abs(fftshift(ifft2(abs(fft2(I, sza, szb)).^2)))./sum(I(:));

x =  (-size(I,2):size(I,2))*ref.PixelExtentInWorldX;
y =  (-size(I,1):size(I,1))*ref.PixelExtentInWorldY;

[X,Y] = ndgrid(x,y);

W = griddedInterpolant(X,Y,C2');
end