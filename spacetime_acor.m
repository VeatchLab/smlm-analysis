function [g,errs,time_edge_cor,N,Norm] = spacetime_acor(x,y,t,tau,r,...
                            smask,tmask,how,timevec,molinframe)
% [G,ERR,TIME_EDGE_COR,N,NORM] = SPACETIME_ACOR(X,Y,T,TAU,R,SMASK,TMASK,...
%                                    HOW,TIMEVEC,MOLINFRAME)
%       space-time autocorrelation function of the points X,Y,T, at TAU and R
%       separations in time and space respectively. HOW specifies whether time-edge-correction
%       should assume uniform density or observed density in time. TIMEVEC and MOLINFRAME
%       are used for the density-corrected (i.e. observed density) time-edge-correction.

    if smask.type ~= 'polygon'
        error('spacetime_acor: mask type not supported')
    end

    maskx = smask.x;
    masky = smask.y;

    T = tmask;

    Dtau = tau(2)-tau(1);
    taubinedges = min(tau)-Dtau/2 : Dtau : max(tau)+Dtau/2;
    Dr = r(2)-r(1);
    rbinedges = min(r)-Dr/2 : Dr : max(r)+Dr/2;
    
    taumin = max(0,min(taubinedges));
    taumax = max(taubinedges);
    rmin = max(0, min(rbinedges));
    rmax = max(rbinedges);
    noutmax = 2e8;
    
    N = closepairs_ts_binned(x,y,t, rmax, numel(r), taumin, taumax, numel(tau));
    
    % basic normalization (no edge corrections)
    area_per_rbin = 2*pi*r'*Dr;
    time_per_tbin = Dtau;
    area = polyarea(maskx, masky);
    duration_excluding_gaps = timewin_duration(T);
    
    density = numel(x)/area/duration_excluding_gaps;
    
    basic_normalization = duration_excluding_gaps*area*density*density*area_per_rbin*time_per_tbin;
    
    edge_cor = spatial_edge_correction(maskx, masky, r);

    if strcmp(how, 'uniform')
        time_edge_cor = time_edge_correction_unif(taubinedges, timevec);
    elseif strcmp(how, 'actual')
       time_edge_cor = time_edge_correction_density(molinframe(:), molinframe(:),...
                                taubinedges, timevec(:), T);
    else
        error('invalid time edge correction method supplied')
    end
    
    g = N./basic_normalization./time_edge_cor./edge_cor;
    
    errs = sqrt(N)./basic_normalization./time_edge_cor./edge_cor;

    Norm = basic_normalization.*time_edge_cor.*edge_cor;
end

