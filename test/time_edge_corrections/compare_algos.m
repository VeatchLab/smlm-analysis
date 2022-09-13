% generate some test data
ll = 1e4;
data = repmat(struct('x', [], 'y', [], 't', []), 1e4, 1);
for i=1:numel(data)
    n = poissrnd(100);
    data(i).x = zeros(1,n);
    data(i).y = zeros(1,n);
    data(i).t = zeros(1,n) + i;
end
x = [data.x];
y = [data.y];
t = [data.t];

timewin = [.5, ll + .5];

taubinedges = [.5:1:ll];

% test that time_edge_correction_density and time_edge_correction_single match
tic;
[tec1] = time_edge_correction_density(t, taubinedges, timewin);
t_density = toc

tic;
[tec2] = time_edge_correction_single(t, taubinedges, timewin);
t_single = toc
