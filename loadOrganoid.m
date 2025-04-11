function [d,si,header] = loadOrganoid(fname,params)


%% Load in, select channel of interest
%%% want the following variables to be constructed: 
%%% d - 
%%% s - 
%%% header - 

[fpath,name,fext] = fileparts(fname);

if strcmpi(fext,'.abf')

    [d,si,header] = abfload(fname); % si is in microseconds
    
    ch_inds = zeros(length(header.recChNames),1);
    for i = 1:length(header.recChNames)
        ch_inds(i) = strcmpi(header.recChNames{i},'Vm1');
    end
    clear i
    assert(sum(ch_inds)==1)
    channel_of_interest = find(ch_inds);
    assert(channel_of_interest<=size(d,2))
    d = d(:,channel_of_interest);

    header.fileType = 'ABF';

elseif strcmpi(fext,'.mat')
    
    g = load(fname);
    [~,varname,~] = fileparts(fname);
    varname = strcat(varname,'_wave_data');
    assert(isfield(g,varname));
    ch_inds = zeros(length(g.(varname).chaninfo),1);
    for ch = 1:length(g.(varname).chaninfo)
        ch_inds(ch) = strcmpi(g.(varname).chaninfo(ch).title,'Vmemb'); % IS THIS ALWAYS TRUE??
    end
    clear ch
    assert(sum(ch_inds)==1)
    channel_of_interest = find(ch_inds);
    
    d = g.(varname).values(:,channel_of_interest);
    si = g.(varname).interval; %This is in seconds!
    si = si*1000000; % convert to microseconds

    g.(varname).values = [];
    header = g.(varname);
    header.fileType = 'Signal-MAT';

%% need to manually multiply by gain as signal outputs something scaled to something else?
    d = d*50;

    clear g varname ch_inds

elseif strcmpi(fext,'.cfs')
    
    fprintf('... File type CFS\n')
    fprintf('... Loading file...')

    % construct params if not specified
    if ~exist('params','var')
        params = struct;
    end
    if ~isfield(params,'concat')
        params.concat = 1; % Concatenate by default
    end
    if ~isfield(params,'showChannels')
    params.showChannels = 0;
    end
    if ~isfield(params,'channelsToLoad')
        params.channelsToLoad = []; % all channels by default
    end
    if ~isfield(params,'manualCheckLoadin')
        params.manualCheckLoadin = 0;
    end

    [d,metadata] = readCFSFile_64b(fname,params);
    
    channels_of_interest = ~isnan(d(1,:,1));

    si = NaN(length(channels_of_interest),1);
    for ch = 1:length(channels_of_interest)
        si(ch) = 1./(metadata.cfsMeta.frameMeta(1).channel(ch).samplingFreq);

    end
    if length(si)>1 % we should sorta check these are the same sampling rates, they really should be...
        for siind = 1:length(si)-2
            assert(si(siind)-si(siind+1)<eps);
        end
    end
    si = si(1)*1000000; %convert to microsec

    metadata.fileType = 'CFS';
    header = metadata;

    clear channels_of_interest metadata

elseif strcmpi(fext,'.edr')

    fprintf('... Data type EDR\n')
    fprintf('... Loading file...\n')

    [data,h] = import_edr(fname);
    si = h.DT;
    si = si*1000000; % convert to microseconds
    ch_inds = zeros(h.NC,1);
    for ch = 1:h.NC
        ch_inds(ch) = strcmpi(h.(strcat('YN',num2str(ch-1))),'Vm1');
    end
    assert(sum(ch_inds)==1)
    channel_of_interest = find(ch_inds)+1; % this +1 is because the first column in d is the timestamps
    
    header = h;
    header.fileType = 'EDR';
    
    fprintf('... Found channel %d, name %s\n',channel_of_interest,h.(strcat('YN',num2str(channel_of_interest-2)))) % -2 because +1 on line 84 and the channels are 0-indexed but matlab is 1-indexed (->+1)
    fprintf('... Unit of recording: %s\n',h.(strcat('YU',num2str(channel_of_interest-2))))

    d = data(:,channel_of_interest);
    
    fprintf('... File loaded successfully!')
    clear data channel_of_interest
else 
    error('No correct file type specified')
end