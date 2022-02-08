function [c, errs, time_edge_cor, N, Norm]= stxcors_from_imagestructs(is, r, tau, how)

% 'actual' means use density-corrected time edge correction
if nargin<4, how = 'actual'; end

nimage = numel(is);
nmask = sum(cellfun(@numel, {is.maskx}));

ii = 0;
% produce a correlation function for each mask in the imagestruct
for i = 1:nimage

    % if data is linked, load it
    if ~isempty(is(i).data)
        data = is(i).data;
    else
        data = load(is(i).data_fname);
    end

    movienum = size(data.data{1}, 1);
    framespermovie = size(data.data{1}, 2);
    nframes = numel(data.data{1});
    
    % put in correct time order
    data.data{1} = reshape(data.data{1}', 1, numel(data.data{1}));
    data.data{2} = reshape(data.data{2}', 1, numel(data.data{2}));
    
    % find the associated 'record', to access relevant metadata
    if isfield(is(i), 'record_fname')
        record_fname = is(i).record_fname;
    else
        [filepath, name, ex] = fileparts(is(i).data_fname);
        record_fname = [filepath '/record.mat'];
    end
    record = load(record_fname);
    
    mdata = record.metadata;
    frame_time = record.metadata.frame_time;
    frames = 1:nframes;
    
    timevec = zeros(nframes, 1);
    moviei_start_time = zeros(movienum, 1);

    % Initialize a timewindow matrix
    T = zeros(0,2);
    for m = 1:movienum
        moviei_start_time(m) = 60 * 60 * 24 * rem(mdata(m).start_time, 1);
        newtimes = frame_time * (1:framespermovie) + ...
            (moviei_start_time(m) - moviei_start_time(1));

        % set the new times in the timevec
        timevec((1:framespermovie) + (m - 1) * framespermovie) = newtimes;

        % add this movie's interval to the time window
        T(m,:) = [newtimes(1) - frame_time/2, newtimes(end) + frame_time/2];
    end

    for j = 1:numel(is(i).maskx)
        ii = ii + 1;

        maskx = is(i).maskx{j};
        masky = is(i).masky{j};
        clear maskeddata
        molinframe1 = zeros(numel(data.data{1}), 1);
        molinframe2 = zeros(numel(data.data{1}), 1);
        for f=1:numel(data.data{1})
            ind = inpolygon(data.data{1}(f).x,data.data{1}(f).y, maskx, masky);
            maskeddata{1}(f).x = data.data{1}(f).x(ind);
            maskeddata{1}(f).y = data.data{1}(f).y(ind);
            maskeddata{1}(f).t = ones(size(data.data{1}(f).x(ind)))*timevec(frames(f));
            molinframe1(f) = sum(ind);
            
            ind = inpolygon(data.data{2}(f).x,data.data{2}(f).y, maskx, masky);
            maskeddata{2}(f).x = data.data{2}(f).x(ind);
            maskeddata{2}(f).y = data.data{2}(f).y(ind);
            maskeddata{2}(f).t = ones(size(data.data{2}(f).x(ind)))*timevec(frames(f));
            molinframe2(f) = sum(ind); 
        end
        
        x1 = [maskeddata{1}(:).x];
        y1 = [maskeddata{1}(:).y];
        x2 = [maskeddata{2}(:).x];
        y2 = [maskeddata{2}(:).y];
        t1 = [maskeddata{1}(:).t];
        t2 = [maskeddata{2}(:).t];
        
        Dtau = tau(2)-tau(1);
        taubinedges = min(tau)-Dtau/2 : Dtau : max(tau)+Dtau/2;
        Dtau/frame_time
        Dr = r(2)-r(1);
        rbinedges = min(r)-Dr/2 : Dr : max(r)+Dr/2;
        
        taumin = min(taubinedges);
        taumax = max(taubinedges);
        rmin = max(0, min(rbinedges));
        rmax = max(rbinedges);
        noutmax = 2e8;
        
        
        %rbinned = crosspairs_rbinned(x1, y1, t1, x2, y2, t2, redges,taumin, taumax, ignore_dr_0)
        
        
        [dx, dy, dt] = crosspairs(x1, y1, t1, x2, y2, t2, rmax, taumin, taumax, noutmax);
        if numel(dx) == noutmax
            warning('Too many pairs! Consider taking a smaller number of data frames or lowering rmax.')
        end
        dr = sqrt(dx.^2+dy.^2);
        N = histcounts2(dr', dt', rbinedges, taubinedges);
        
        % basic normalization (no edge corrections)
        area_per_rbin = 2*pi*r'*Dr;
        time_per_tbin = Dtau;
        area = polyarea(maskx, masky);
        total_duration = timevec(frames(end))-timevec(frames(1));
        duration_excluding_gaps = timewin_duration(T);
        
        density1 = numel(x1)/area/duration_excluding_gaps;
        density2 = numel(x2)/area/duration_excluding_gaps;
        
        basic_normalization = duration_excluding_gaps*area*density1*density2*area_per_rbin*time_per_tbin;
        
        edge_cor = edge_correction(maskx, masky, r);
        
        if strcmp(how, 'uniform')
            time_edge_cor = time_edge_correction_unif(taubinedges, timevec);
        elseif strcmp(how, 'actual')
            time_edge_cor = time_edge_correction(molinframe1, molinframe2,...
                                taubinedges, timevec, T);
        else
            error('invalid time edge correction method supplied')
        end
        % kluge until this is fixed in the time_edge_correction()s.
        %time_edge_cor = time_edge_cor/(Dtau/frame_time);
        
        c(ii, :, :) = N./basic_normalization./time_edge_cor./edge_cor;
        
        errs(ii, :, :) = sqrt(N)./basic_normalization./time_edge_cor./edge_cor;
        
        Norm = basic_normalization.*time_edge_cor.*edge_cor;
    end
end
end


function correction = edge_correction(maskx, masky, r)
% isotropic spatial edge correction
area = polyarea(maskx, masky);

% compute bin geometry
binwidth = diff([0 r]);
rcenter = r - binwidth/2; % find bin centers

% Compute an edge-correction factor, averaged over thetas
ntheta = 30;
thetas = pi*(1:ntheta)/ntheta;
correction = zeros(size(r));
for i = 1:numel(r)
    for j = 1:ntheta
        dx = cos(thetas(j))*rcenter(i);
        dy = sin(thetas(j))*rcenter(i);
        correction(i) = correction(i) + wij(maskx, masky, dx, dy)/ntheta;
    end
end
correction = correction'/area;
end

function w = wij(maskx, masky, dx, dy)
[x, y] = polybool('intersection', maskx, masky, maskx + dx, masky + dy);
w = polyarea(x,y);
end

function taufactor = time_edge_correction(Nperframe1, Nperframe2, tau, ...
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

timediffs = timevec - timevec';
weights = Nperframe1.*Nperframe2';% is this direction right?
% count = 1;
% for i = 1:tmax
%     for j = 1:tmax
%         timediffs(count) = timevec(j) - timevec(i);
%         weights(count) = Nperframe1(i)*Nperframe2(j);
%         count = count + 1;
%     end
% end

dtau = diff(tau);

[~, ~, bin] = histcounts(timediffs, tau);
inds = bin>0; % & bin <= (tau(end)-tau(1))/dt;
exptauperbin = accumarray(bin(inds), weights(inds)')./dtau(:);
taufactor = exptauperbin*timewin_duration(timewin)/(sum(Nperframe1)*sum(Nperframe2));
taufactor = taufactor';

end

function taufactor = time_edge_correction_unif(tau, timevec)
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
