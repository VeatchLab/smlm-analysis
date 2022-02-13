function taufactor = time_edge_correction_unif(tau, timewin)
tmax = numel(timevec);
timediffs = zeros(tmax*(tmax-1)/2, 1);
count = 1;
for i = 1:tmax-1
    for j = i:tmax
        timediffs(count) = timevec(j) - timevec(i);
        count = count + 1;
    end
end

dt = tau(2)-tau(1);

[~, ~, bin] = histcounts(timediffs, [tau]);% tau(end)+2*dt]);
inds = bin>0 & bin <= tau(end)/dt + 1;
exptauperbin = accumarray(bin(inds), 1);
taufactor = exptauperbin/numel(timevec);
taufactor = taufactor';

end

