function l = timewin_duration(T)
% L = TIMEWIN_DURATION(T) find the duration of a temporal window
%       T is an ordered, disjoint set of N closed intervals,
%       represented as an Nx2 matrix, so that T(i,1) and T(i,2)
%       are the start and end times of the ith interval, respectively

% check that T is a valid time window
if ~timewin_isvalid(T)
    error('timewin_duration: an invalid time window T was provided');
end

% differences along second dimension, so (end - start) of each interval
% then just add 'em up.
l = sum(diff(T,1,2));
