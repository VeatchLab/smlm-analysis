function [c, errs, time_edge_cor, N, Norm]= stxcors_from_imagestructs(is, r, tau, how)

if nargin<4, how = 'actual'; end



nimage = numel(is);
nmask = sum(cellfun(@numel, {is.maskx}));
%c = zeros(nmask, numel(r));
%errs = zeros(nmask, numel(r));


ii = 0;
for i = 1:nimage

    if ~isempty(is(i).data)
        data = is(i).data;
    else
        data = load(is(i).data_fname);
    end
    
    framespermovie = size(data.data{1}, 2);
    
    % put in correct time order
    data.data{1} = reshape(data.data{1}', 1, numel(data.data{1}));
    data.data{2} = reshape(data.data{2}', 1, numel(data.data{2}));
    
    if isfield(is(i), 'record_fname')
        record_fname = is(i).record_fname;
    else
        [filepath, name, ex] = fileparts(is(i).data_fname);
        record_fname = [filepath '/record.mat'];
    end
    record = load(record_fname);
    
    mdata = record.metadata;
    frame_time = record.metadata.frame_time;
    nframes = numel(data.data{1});
    frames = 1:nframes;
    
    %framespermovie = size(data.data{1}, 2);
    movienum = numel(data.data{1})/framespermovie;
    timevec = zeros(nframes, 1);
    moviei_start_time = zeros(movienum, 1);
    for m = 1:movienum
        moviei_start_time(m) = 60 * 60 * 24 * rem(mdata(m).start_time, 1);
        timevec((1:framespermovie) + (m - 1) * framespermovie) = frame_time * (1:framespermovie) + ...
            (moviei_start_time(m) - moviei_start_time(1));
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
        duration_excluding_gaps = frame_time*numel(timevec);
        
        density1 = numel(x1)/area/duration_excluding_gaps;
        density2 = numel(x2)/area/duration_excluding_gaps;
        
        basic_normalization = duration_excluding_gaps*area*density1*density2*area_per_rbin*time_per_tbin;
        
        edge_cor = edge_correction(maskx, masky, r);
        %time_edge_cor = time_edge_correction(molinframe1, molinframe2, taubinedges, timevec);
        %time_edge_cor = (duration-tau)/duration;
        
        if strcmp(how, 'uniform')
            time_edge_cor = time_edge_correction_unif(taubinedges, timevec);
        elseif strcmp(how, 'actual')
            time_edge_cor = time_edge_correction(molinframe1, molinframe2, taubinedges, timevec);
        else
            error('invalad time edge correction method supplied')
        end
	time_edge_cor = time_edge_cor/(Dtau/frame_time);
        
        
        
        c(ii, :, :) = N./basic_normalization./time_edge_cor./edge_cor;
        
        
        errs(ii, :, :) = sqrt(N)./basic_normalization./time_edge_cor./edge_cor;
        
        Norm = basic_normalization.*time_edge_cor.*edge_cor;
    end
end
end


function correction = edge_correction(maskx, masky, r)

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

function taufactor = time_edge_correction(Nperframe1, Nperframe2, tau, timevec)

%tmax = numel(Nperframe1);

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

dt = tau(2)-tau(1);

[~, ~, bin] = histcounts(timediffs, tau);% tau(end)+2*dt]);
inds = bin>0 & bin <= (tau(end)-tau(1))/dt;
exptauperbin = accumarray(bin(inds), weights(inds)');
taufactor = exptauperbin/mean(Nperframe1)/mean(Nperframe2)/numel(timevec);
taufactor = taufactor';

end




% 
% 
% function taufactor = time_edge_correction(Nperframe1, Nperframe2, tau, timevec)
% 
% tmax = numel(Nperframe1);
% 
% timediffs = zeros(tmax*(tmax-1)/2, 1);
% weights = zeros(tmax*(tmax-1)/2, 1);
% count = 1;
% for i = 1:tmax-1
%     for j = i:tmax
%         timediffs(count) = timevec(j) - timevec(i);
%         weights(count) = Nperframe1(i)*Nperframe2(j);
%         count = count + 1;
%     end
% end
% 
% dt = tau(2)-tau(1);
% 
% [~, ~, bin] = histcounts(timediffs, [tau]);% tau(end)+dt]);
% inds = bin>0 & bin <= tau(end)/dt + 1;
% exptauperbin = accumarray(bin(inds), weights(inds)');
% taufactor = exptauperbin/mean(Nperframe1)/mean(Nperframe2)/numel(timevec);
% taufactor = taufactor';
% 
% end

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


