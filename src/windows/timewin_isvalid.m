function valid = timewin_isvalid(T)
% L = TIMEWIN_ISVALID(T) check that the time window is valid/allowable
%       T should be an ordered, disjoint set of N closed intervals,
%       represented as an Nx2 matrix, so that T(i,1) and T(i,2)
%       are the start and end times of the ith interval, respectively

goodsz = size(T,2) == 2;
long = T';
long = long(:);
% comparison should be strict (i.e. not allowing equal values)
goodsort = issorted(long,'strictascend');

valid = goodsz & goodsort;
