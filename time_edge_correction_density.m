function taufactor = time_edge_correction_density(Nperframe1, Nperframe2, tau, ...
                            timevec, timewin)
% density corrected time edge correction
tmax = numel(Nperframe1);

% compare all sizes to the size of timevec to make sure they match
tsz = size(timevec);
if ~isequal(tsz, size(Nperframe1)) || ~isequal(tsz, size(Nperframe2))
    error('time_edge_correction: size of Nperframe''s must match timevec');
end

% check validity of timewin argument
if ~timewin_isvalid(timewin)
    error('time_edge_correction: invalid time window provided');
end

timediffs = zeros(tmax*(tmax-1)/2, 1);
weights = zeros(tmax*(tmax-1)/2, 1);
count = 1;
for i = 1:tmax-1
    for j = i:tmax
        timediffs(count) = timevec(j) - timevec(i);
        weights(count) = Nperframe1(i)*Nperframe2(j);
        count = count + 1;
    end
end

dtau = diff(tau);

[~, ~, bin] = histcounts(timediffs, tau);
inds = bin>0;
exptauperbin = accumarray(bin(inds), weights(inds)')./dtau(:);
taufactor = exptauperbin*timewin_duration(timewin)/(sum(Nperframe1)*sum(Nperframe2));
taufactor = taufactor';

end
