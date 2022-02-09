function correction = spatial_edge_correction(maskx, masky, r)
% isotropic spatial edge correction
area = polyarea(maskx, masky);

% compute bin geometry
binwidth = diff([0 r]);
rcenter = r - binwidth/2; % find bin centers

% Compute an edge-correction factor, averaged over thetas
ntheta = 30;
thetas = pi*(1:ntheta)/ntheta;
correction = zeros(size(r));
for i = 1:numel(r)
    for j = 1:ntheta
        dx = cos(thetas(j))*rcenter(i);
        dy = sin(thetas(j))*rcenter(i);
        correction(i) = correction(i) + wij(maskx, masky, dx, dy)/ntheta;
    end
end
correction = correction'/area;
end

function w = wij(maskx, masky, dx, dy)
[x, y] = polybool('intersection', maskx, masky, maskx + dx, masky + dy);
w = polyarea(x,y);
end

