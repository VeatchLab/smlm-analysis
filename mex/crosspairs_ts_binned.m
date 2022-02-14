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
