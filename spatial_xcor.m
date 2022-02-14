function [g,errs,N,Norm] = spatial_xcor(x1,y1,x2,y2,r,smask)
% [G,ERR,N,NORM] = SPATIAL_XCOR(X,Y,R,SMASK)
%       spatial autocorrelation function of the points X,Y, at R
%       separations in space.

    if smask.type ~= 'polygon'
        error('spacetime_acor: mask type not supported')
    end

    maskx = smask.x;
    masky = smask.y;

    Dr = r(2)-r(1);
    rbinedges = min(r)-Dr/2 : Dr : max(r)+Dr/2;
    
    rmin = max(0, min(rbinedges));
    rmax = max(rbinedges);

    % to pass to closepairs
    t1 = zeros(size(x1));
    t2 = zeros(size(x2));
    taumin = -1;
    taumax = 1;
    
    N = crosspairs_ts_binned(x1,y1,t1,x2,y2,t2,rmax, numel(r), taumin, taumax, 1);
    
    % basic normalization (no edge corrections)
    area_per_rbin = 2*pi*r'*Dr;
    area = polyarea(maskx, masky);
    
    density1 = numel(x1)/area;
    density2 = numel(x2)/area;
    basic_normalization = area*density1*density2*area_per_rbin;
    edge_cor = spatial_edge_correction(maskx, masky, r);

    g = N./basic_normalization./edge_cor;
    
    errs = sqrt(N)./basic_normalization./edge_cor;

    Norm = basic_normalization.*edge_cor;
end
