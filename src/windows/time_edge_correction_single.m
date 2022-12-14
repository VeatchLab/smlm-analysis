function [taufactor,norm] = time_edge_correction_single(t,tau_edges,timewin)
% TIME_EDGE_CORRECTION_SINGLE temporal edge correction for correlation
% function

% Copyright (C) 2021 Thomas Shaw, and Sarah Veatch
% This file is part of SMLM SPACETIME RESOLUTION
% SMLM SPACETIME RESOLUTION is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% SMLM SPACETIME RESOLUTION is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with SMLM SPACETIME RESOLUTION.  If not, see <https://www.gnu.org/licenses/>

ntout = numel(tau_edges) -1;
% check validity of timewin argument
if ~timewin_isvalid(timewin) || size(timewin,1) > 1
    error('time_edge_correction_single: invalid time window provided');
end

ts = sort(unique(t));
dt = min(diff(ts));
tbinedges = (min(ts) - dt/2):dt:(max(ts) + dt);
%fprintf('There are %f times as many points as unique times\n', numel(t)/numel(timevec));

Nperframe = histcounts(t, tbinedges);

weights = conv(Nperframe, fliplr(Nperframe));
nt = round((max(ts) - min(ts))/dt);
timediffs = (-nt:nt)*dt;

[~, ~, bin] = histcounts(timediffs, tau_edges);
inds = bin>0;
exptauperbin = accumarray(bin(inds)', weights(inds),[ntout,1])/dt;

norm = numel(t)^2 / timewin_duration(timewin);
taufactor = exptauperbin'/norm;
