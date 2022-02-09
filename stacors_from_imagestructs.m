function [c, errs, time_edge_cor, N, Norm] = stacors_from_imagestructs(is, r, tau, how)

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

    % find the associated 'record', to access relevant metadata
    if isfield(is(i), 'record_fname')
        record_fname = is(i).record_fname;
    else
        [filepath, name, ex] = fileparts(is(i).data_fname);
        record_fname = [filepath '/record.mat'];
    end

    if strcmp(record_fname, 'UseDefaultTiming')
        % 'UseDefaultTiming' means frames are all sequential and separated by 1
        frame_time = 1;
        frames = 1:nframes;
        timevec = frames;
        % timewindow matrix is just one interval
        T = [timevec(1) - frame_time/2, timevec(end) + frame_time/2];
    else
        % Otherwise assemble timing data from the metadata that is
        % loaded into record
        record = load(record_fname, 'metadata');
        
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
    end
    
    % do autocorrelation for each channel present in the data
    for channel=1:numel(data.data)
        % put in correct time order
        d = reshape(data.data{channel}', 1, numel(data.data{channel}));
        
        for j = 1:numel(is(i).maskx)
            ii = ii + 1;
            
            maskx = is(i).maskx{j};
            masky = is(i).masky{j};

            for f=1:numel(d)
                ind = inpolygon(d(f).x,d(f).y, maskx, masky);
                maskeddata(f).x = d(f).x(ind);
                maskeddata(f).y = d(f).y(ind);
                maskeddata(f).t = ones(size(d(f).x(ind)))*timevec(frames(f));
            end

            molinframe = arrayfun(@(s) numel(s.x), maskeddata);
            
            x = [maskeddata(:).x];
            y = [maskeddata(:).y];
            t = [maskeddata(:).t];
            
            % call out to autocorrelation function calculation
            smask = struct('x', maskx, 'y', masky, 'type', 'polygon');
            [cc,ee,time_edge_cor, N, Norm] = spacetime_acor(x,y,t,...
                                        tau,r,smask, T, how, timevec, molinframe);

            c{channel}(i,:,:) = cc;
            errs{channel}(i,:,:) = ee;
        end
    end
end
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
