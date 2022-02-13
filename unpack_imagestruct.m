function out = unpack_imagestruct(is)

nimage = numel(is);

% Handle case of multiple images
% if there's more than one image, handle each one separately
if nimage > 1
    out = unpack_imagestruct(is(1));
    [out(:).imageid] = deal(1);
    for i=2:nimage
        % each image may produce multiple outputs
        iout = unpack_imagestruct(is(i));
        [iout(:).imageid] = deal(i);
        out = [out, iout];
    end
    return
end

ii = 0;

% if data is linked, load it
if ~isempty(is.data)
    data = is.data;
else
    data = load(is.data_fname);
end

% figure out how long movies are, and how many
movienum = size(data.data{1}, 1);
framespermovie = size(data.data{1}, 2);
nframes = numel(data.data{1});

% find the associated 'record', to access relevant metadata
if isfield(is, 'record_fname')
    record_fname = is.record_fname;
else
    [filepath, name, ex] = fileparts(is.data_fname);
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
    
% do unpacking for each channel present in the data
for channel=1:numel(data.data)
    % put in correct time order
    d = reshape(data.data{channel}', 1, numel(data.data{channel}));
    
    for j = 1:numel(is.maskx)
        ii = ii + 1;
        
        maskx = is.maskx{j};
        masky = is.masky{j};

        for f=1:numel(d)
            ind = inpolygon(d(f).x,d(f).y, maskx, masky);
            maskeddata(f).x = d(f).x(ind);
            maskeddata(f).y = d(f).y(ind);
            maskeddata(f).t = ones(size(d(f).x(ind)))*timevec(frames(f));
        end

        x = [maskeddata(:).x];
        y = [maskeddata(:).y];
        t = [maskeddata(:).t];
        
        % Make struct for spatial mask
        smask = struct('x', maskx, 'y', masky, 'type', 'polygon');

        % desired fields:
        % x,y,t,spacewin,timewin,frame_time,timevec,channel,maskid
        out(ii) = struct('x', x, 'y', y, 't', t, 'spacewin', smask, 'timewin', T,...
            'channel', channel, 'maskid', j,...
            'imageid', 1); % overwritten later if there is more than one
    end
end
end