function [c, errs, time_edge_cor, N, Norm] = stacors_from_imagestructs(is, r, tau, how)

if nargin<4, how = 'actual'; end


nimage = numel(is);
nmask = sum(cellfun(@numel, {is.maskx}));
%c = zeros(nmask, numel(r),numel(r));
%errs = zeros(nmask, numel(r),numel(r));


ii = 0;
for i = 1:nimage

    if ~isempty(is(i).data)
        data = is(i).data;
    else
        data = load(is(i).data_fname);
    end
    
    framespermovie = size(data.data{1}, 2);
    if isfield(is(i), 'record_fname')
        record_fname = is(i).record_fname;
    else
        [filepath, name, ex] = fileparts(is(i).data_fname);
        record_fname = [filepath '/record.mat'];
    end
    
    if ~strcmp(record_fname, 'UseDefaultTiming')
        record = load(record_fname, 'metadata');
        
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
        
    else
        frame_time = 1;
        nframes = numel(data.data{1});
        frames = 1:nframes;
        timevec = frames;
    end
    
    for channel=1:numel(data.data)
        
        % put in correct time order
        d = reshape(data.data{channel}', 1, numel(data.data{channel}));
        
        
        for j = 1:numel(is(i).maskx)
            ii = ii + 1;
            
            maskx = is(i).maskx{j};
            masky = is(i).masky{j};
            clear maskeddata
            molinframe = zeros(numel(data), 1);
            for f=1:numel(d)
                ind = inpolygon(d(f).x,d(f).y, maskx, masky);
                maskeddata(f).x = d(f).x(ind);
                maskeddata(f).y = d(f).y(ind);
                maskeddata(f).t = ones(size(d(f).x(ind)))*timevec(frames(f));
                molinframe(f) = sum(ind);
            end
            
            x = [maskeddata(:).x];
            y = [maskeddata(:).y];
            t = [maskeddata(:).t];
            
            Dtau = tau(2)-tau(1);
            taubinedges = min(tau)-Dtau/2 : Dtau : max(tau)+Dtau/2;
            Dr = r(2)-r(1);
            rbinedges = min(r)-Dr/2 : Dr : max(r)+Dr/2;
            
            taumin = max(0,min(taubinedges));
            taumax = max(taubinedges);
            rmin = max(0, min(rbinedges));
            rmax = max(rbinedges);
            noutmax = 2e8;
            
            N = closepairs_ts_binned(x,y,t, rmax, numel(r), taumin, taumax, numel(tau));
            
            % basic normalization (no edge corrections)
            area_per_rbin = 2*pi*r'*Dr;
            time_per_tbin = Dtau;
            area = polyarea(maskx, masky);
            %total_duration = timevec(frames(end))-timevec(frames(1));
            duration_excluding_gaps = frame_time*numel(timevec);
            
            density = numel(x)/area/duration_excluding_gaps;
            
            basic_normalization = duration_excluding_gaps*area*density*density*area_per_rbin*time_per_tbin;
            
            edge_cor = edge_correction(maskx, masky, r);
            %time_edge_cor = time_edge_correction(molinframe, molinframe, taubinedges, timevec);

            if strcmp(how, 'uniform')
                time_edge_cor = time_edge_correction_unif(taubinedges, timevec);
            elseif strcmp(how, 'actual')
               time_edge_cor = time_edge_correction(molinframe, molinframe, taubinedges, timevec);
            else
                error('invalad time edge correction method supplied')
            end
                
                
            %time_edge_cor = (duration-tau)/duration;
            %time_edge_cor
            
            cc = N./basic_normalization./time_edge_cor./edge_cor;
            c{channel}(i, :, :) = cc;
            
            errs{channel}(i, :, :) = sqrt(N)./basic_normalization./time_edge_cor./edge_cor;

	    Norm = basic_normalization.*time_edge_cor.*edge_cor;
            
        end
        
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

tmax = numel(Nperframe1);

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

dt = tau(2)-tau(1);

[~, ~, bin] = histcounts(timediffs, [tau]);% tau(end)+2*dt]);
inds = bin>0 & bin <= tau(end)/dt + 1;
exptauperbin = accumarray(bin(inds), weights(inds)');
taufactor = exptauperbin/mean(Nperframe1)/mean(Nperframe2)/numel(timevec);
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
