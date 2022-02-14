function [counts] = closepairs_ts_binned(x, y, t, rmax, nrout, taumin, taumax, ntout)

% Sort based on x
[x, s] = sort(x(:));

t = t(s);
y = y(s);

[counts] = Fclosepairs_ts_binned(x,y,t, rmax, uint64(nrout), taumin, taumax, uint64(ntout));

counts = double(counts);
