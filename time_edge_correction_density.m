function [taufactor,norm] = time_edge_correction_density(t,tau_edges,timewin)
ntout = numel(tau_edges) -1;
% check validity of timewin argument
if ~timewin_isvalid(timewin)
    error('time_edge_correction: invalid time window provided');
end

timevec = sort(unique(t));
%fprintf('There are %f times as many points as unique times\n', numel(t)/numel(timevec));

Nperframe = arrayfun(@(tt) sum(t == tt), timevec);

tmax = numel(timevec);

timediffs = zeros(tmax*(tmax-1)/2, 1);
weights = zeros(tmax*(tmax-1)/2, 1);
count = 1;
for i = 1:tmax-1
    for j = i:tmax
        timediffs(count) = timevec(j) - timevec(i);
        weights(count) = Nperframe(i)*Nperframe(j);
        count = count + 1;
    end
end

dtau = diff(tau_edges);

[~, ~, bin] = histcounts(timediffs, tau_edges);
inds = bin>0;
exptauperbin = accumarray(bin(inds), weights(inds)',[ntout,1])./dtau(:);

norm = numel(t)^2 / timewin_duration(timewin);
taufactor = exptauperbin'/norm;
