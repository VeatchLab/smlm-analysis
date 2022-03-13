% Copyright (C) 2021 Thomas Shaw, Frank Fazekas and Sarah Veatch
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

function [counts] = crosspairs_ts_binned(x1, y1, t1, x2, y2, t2, ...
        rmax, nrout, taumin, taumax, ntout)

% Sort based on x
[x1, s1] = sort(x1(:));
[x2, s2] = sort(x2(:));

t1 = t1(s1);
y1 = y1(s1);
t2 = t2(s2);
y2 = y2(s2);

[counts] = Fcrosspairs_ts_binned(x1,y1,t1, x2, y2, t2,...
                                    rmax, uint64(nrout), taumin, taumax, uint64(ntout));

counts = double(counts);
