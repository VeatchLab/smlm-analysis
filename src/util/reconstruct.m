function I = reconstruct(data, iref)
% use iref

left = iref.XWorldLimits(1);
right = iref.XWorldLimits(2);
top = iref.YWorldLimits(1);
bottom = iref.YWorldLimits(2);
pwidth = iref.PixelExtentInWorldX;
pheight = iref.PixelExtentInWorldY;

xedges = left:pwidth:right;
yedges = top:pheight:bottom;

xs = [data.x]; ys = [data.y];

I = histcounts2(ys, xs, yedges, xedges);
