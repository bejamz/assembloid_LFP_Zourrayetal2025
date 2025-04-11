function [] = analyseOrganoidLFP(xlsfile,main_params)


%% Check inputs and parse params

assert(exist(xlsfile,'file')==2) % need the file to exist!

if ~exist('main_params','var')
    main_params = struct;
end
if ~isfield(main_params,'desired_rate')
    main_params.desired_rate = 600;
end
if ~isfield(main_params,'block_length_sec')
    main_params.block_length_sec = 15;
    fprintf('Setting default block length: %d seconds \n',main_params.block_length_sec);
end
if ~isfield(main_params,'block_length_min')
    main_params.block_length_min = 5;
end

segment_length_min = 30; % NEED TO CHANGE LOGIC BELOW BEFORE CHANGING THIS VALUE

%% Define frequency bands of interest

if ~isfield(main_params,'deltaBand')
    main_params.run_delta = 1;
    main_params.deltaBand = [1,4];
end
if ~isfield(main_params,'thetaBand')
    main_params.run_theta = 1;
    main_params.thetaBand = [5,12];
end
if ~isfield(main_params,'betaBand')
    main_params.run_beta = 1;
    main_params.betaBand = [15,25];
end
if ~isfield(main_params,'gammaBand')
    main_params.run_gamma = 1;
    main_params.gammaBand = [25,100];
end
if ~isfield(main_params,'slowgammaBand')
    main_params.run_sgamma = 1;
    main_params.slowgammaBand = [25,50];
end
if ~isfield(main_params,'fastgammaBand')
    main_params.run_fgamma = 1;
    main_params.fastgammaBand = [50,100];
end
if ~isfield(main_params,'rippleBand')
    main_params.run_ripples = 1;
    main_params.rippleBand = [150,250];
end

params = struct;

%% Read in xls file

[~,~,raw] = xlsread(xlsfile);

%% Preallocate arrays for analysis at end

group = NaN(size(raw,1),1);

allpower = NaN(size(raw,1),4);
thetapower = NaN(size(raw,1),4);
sgammapower = NaN(size(raw,1),4);
fgammapower = NaN(size(raw,1),4);

%% Extract cols from xls file

for col = 1:size(raw,2)
    switch raw{2,col}
        case 'orgid'
            orgCol = col;
        case 'group' % Group is for e.g. control or patient
            groupCol = col;
        case 'protocol'
            protocolCol = col;
        case 'protocol_notes'
            protocolNotesCol = col;
        case 'filepath'
            fpathCol = col;
        case 'filename'
            fnameCol = col;
        case 'kcl_timestamp_type'
            kclTsTypeCol = col;
        case 'kcl_time_sec'
            kclTimeCol = col; % = baseline time, if there's a perturbation!
        case 'kcl_marker'
            kclMarkerCol = col;
        case 'include'
            includeCol = col;
        case 'filetype'
            filetypeCol = col;
        case 'saveFile'
            saveFileCol = col;
        case 'line'
            differentiationCol = col;
        case 'date'
            dateCol = col;
        case 'batch' 
            batchCol = col;
        case 'age'
            ageCol = col;
    end
end



%% Main for loop
for org = 3:size(raw,1)

    if isnan(raw{org,orgCol})
        continue;
    end

    if ~raw{org,includeCol} % if not due to be included in the dataset
        continue;
    end

    fpath = raw{org,fpathCol};
    fname = raw{org,fnameCol};
    fname = [fpath,'\',fname];

    assert(exist(fname,'file')==2)
    
    group(org) = raw{org,groupCol};

    organoid = struct;
    % put metadata into the struct
    organoid.Meta.OrgName = raw{org,orgCol};
    organoid.Meta.Group = raw{org,groupCol};

    [data,si,meta] = loadOrganoid(fname,params);
    data = data(:,1); % remove the current one to halve data use.  This should be changed a bit!

    [data,fs] = downsampleData(data,si,main_params.desired_rate); % data is often obscenely huge, so let's overwrite it here for memory reasons
    si = 1000000/fs; %in microseconds

    [data] = preprocessLFP(data,si,params);

    organoid.Data.data = data;
    organoid.Data.Fs = fs;
    organoid.Data.si = si;
    organoid.Data.meta = meta;

    %% Run Hilbert analyses

    bandpassTheta = designfilt('bandpassiir','FilterOrder',20, ...
             'HalfPowerFrequency1',main_params.thetaBand(1),'HalfPowerFrequency2',main_params.thetaBand(2), ...
             'SampleRate',fs);

    bandpassGamma = designfilt('bandpassiir','FilterOrder',20, ...
             'HalfPowerFrequency1',main_params.gammaBand(1),'HalfPowerFrequency2',main_params.gammaBand(2), ...
             'SampleRate',fs);

    bandpassDelta = designfilt('bandpassiir','FilterOrder',20, ...
             'HalfPowerFrequency1',main_params.deltaBand(1),'HalfPowerFrequency2',main_params.deltaBand(2), ...
             'SampleRate',fs);

    bandpassRipple = designfilt('bandpassiir','FilterOrder',20, ...
             'HalfPowerFrequency1',main_params.rippleBand(1),'HalfPowerFrequency2',main_params.rippleBand(2), ...
             'SampleRate',fs);

    d_delta = filtfilt(bandpassDelta,data(:,1));
    d_theta = filtfilt(bandpassTheta,data(:,1));
    d_gamma = filtfilt(bandpassGamma,data(:,1));
    d_ripple = filtfilt(bandpassRipple,data(:,1));

    [inst_amp_delta,inst_phase_delta] = processHilberts(d_delta,si);
    [inst_amp_theta,inst_phase_theta] = processHilberts(d_theta,si);
    [inst_amp_gamma,inst_phase_gamma] = processHilberts(d_gamma,si);

    [phaseamps,phase_bins] = phaseAmplitudeCoupling(inst_phase_theta,inst_amp_gamma);

    organoid.HilbertAnalysis.inst_amp_delta = inst_amp_delta;
    organoid.HilbertAnalysis.inst_phase_delta = inst_phase_delta;
    organoid.HilbertAnalysis.inst_amp_theta = inst_amp_theta;
    organoid.HilbertAnalysis.inst_phase_theta = inst_phase_theta;
    organoid.HilbertAnalysis.inst_amp_gamma = inst_amp_gamma;
    organoid.HilbertAnalysis.inst_phase_gamma = inst_phase_gamma;
    organoid.HilbertAnalysis.ThetaGammaPhaseAmps = phaseamps;
    organoid.HilbertAnalysis.ThetaGammaPhaseBins = phase_bins;
    
    [phaseamps,phase_bins] = phaseAmplitudeCoupling(inst_phase_delta,inst_amp_gamma);
    organoid.HilbertAnalysis.DeltaGammaPhaseAmps = phaseamps;
    organoid.HilbertAnalysis.DeltaGammaPhaseBins = phase_bins;
    if exist('circ_r')==2
        organoid.HilbertAnalysis.DeltaGammaR = circ_r(deg2rad(phase_bins(:)),phaseamps(:));
    end
    
    clear phaseamps phase_bins inst_amp_delta inst_phase_delta inst_amp_theta inst_phase_theta inst_amp_gamma inst_phase_gamma d_delta d_theta d_gamma d_ripple

    %% Short block analysis - %d SECONDS

    % Standard statistics and spectral measurements
    [stds_smallblocks,norm_coastline_smallblocks,coastline_smallblocks,rms_smallblocks,spectralstuff] = shortBlockAnalysis(data(:,1),fs,main_params.block_length_sec,params);
    
    % Bandpowers
    [short_bandpow_all,~,~]                             = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_sec);
    [~,short_bandpow_delta,short_proportion_delta]      = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_sec/60,main_params.deltaBand,params);
    [~,short_bandpow_theta,short_proportion_theta]      = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_sec/60,main_params.thetaBand,params);
    [~,short_bandpow_beta,short_proportion_beta]        = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_sec/60,main_params.betaBand,params);
    [~,short_bandpow_sgamma,short_proportion_sgamma]    = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_sec/60,main_params.slowgammaBand,params);
    [~,short_bandpow_fgamma,short_proportion_fgamma]    = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_sec/60,main_params.fastgammaBand,params);
    [~,short_bandpow_hfos,short_proportion_hfos]        = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_sec/60,main_params.rippleBand,params);
    
    organoid.SmallBlockAnalysis.Meta.SmallBlockLengthSec    = main_params.block_length_sec;
    organoid.SmallBlockAnalysis.stdevs_smallblocks          = stds_smallblocks; % DO NOT DELETE - THIS IS USED LATER
    organoid.SmallBlockAnalysis.norm_coastline_smallblocks  = norm_coastline_smallblocks;
    organoid.SmallBlockAnalysis.coastline_smallblocks       = coastline_smallblocks;
    organoid.SmallBlockAnalysis.rms_smallblocks             = rms_smallblocks;
    organoid.SmallBlockAnalysis.bandpow_all                 = short_bandpow_all;
    organoid.SmallBlockAnalysis.bandpow_delta               = short_bandpow_delta;
    organoid.SmallBlockAnalysis.bandpow_theta               = short_bandpow_theta;
    organoid.SmallBlockAnalysis.bandpow_beta                = short_bandpow_beta;
    organoid.SmallBlockAnalysis.bandpow_sgamma              = short_bandpow_sgamma;
    organoid.SmallBlockAnalysis.bandpow_fgamma              = short_bandpow_fgamma;
    organoid.SmallBlockAnalysis.bandpow_hfos                = short_bandpow_hfos;
    organoid.SmallBlockAnalysis.proportion_delta            = short_proportion_delta;
    organoid.SmallBlockAnalysis.proportion_theta            = short_proportion_theta;
    organoid.SmallBlockAnalysis.proportion_beta             = short_proportion_beta;
    organoid.SmallBlockAnalysis.proportion_sgamma           = short_proportion_sgamma;
    organoid.SmallBlockAnalysis.proportion_fgamma           = short_proportion_fgamma;
    organoid.SmallBlockAnalysis.proportion_hfos             = short_proportion_hfos;
    organoid.SmallBlockAnalysis.SpectralMeasurements        = spectralstuff;
    
    clear short_bandpow_all short_bandpow_delta short_bandpow_theta short_bandpow_beta short_bandpow_sgamma short_bandpow_fgamma short_bandpow_hfos ...
        short_proportion_delta short_proportion_theta short_proportion_beta short_proportion_sgamma short_proportion_fgamma short_proportion_hfos
    
    %% Large blocks - %d MINUTES
    
    % Standard statistics and spectral measurements
    [stds_largeblocks,norm_coastline_largeblocks,coastline_largeblocks] = shortBlockAnalysis(data(:,1),fs,main_params.block_length_min,params);
    
    % Bandpowers
    [bandpow_all,~,~]                       = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_min);
    [~,bandpow_delta,proportion_delta]      = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_min,main_params.deltaBand,params);
    [~,bandpow_theta,proportion_theta]      = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_min,main_params.thetaBand,params);
    [~,bandpow_sgamma,proportion_sgamma]    = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_min,main_params.slowgammaBand,params);
    [~,bandpow_fgamma,proportion_fgamma]    = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_min,main_params.fastgammaBand,params);
    [~,bandpow_hfos,proportion_hfos]        = bandPowerAnalysisBlocks(data(:,1),fs,main_params.block_length_min,main_params.rippleBand,params);

    organoid.LargeBlockAnalysis.Meta.LargeBlockLengthMin        = main_params.block_length_min;
    organoid.LargeBlockAnalysis.stds_largeblocks                = stds_largeblocks;
    organoid.LargeBlockAnalysis.coastline_largeblocks           = coastline_largeblocks;
    organoid.LargeBlockAnalysis.norm_coastline_largeblocks      = norm_coastline_largeblocks;
    organoid.LargeBlockAnalysis.BandpowerAll                    = bandpow_all;
    organoid.LargeBlockAnalysis.bandpow_delta                   = bandpow_delta;
    organoid.LargeBlockAnalysis.bandpow_theta                   = bandpow_theta;
    organoid.LargeBlockAnalysis.bandpow_sgamma                  = bandpow_sgamma;
    organoid.LargeBlockAnalysis.bandpow_fgamma                  = bandpow_fgamma;
    organoid.LargeBlockAnalysis.bandpow_hfos                    = bandpow_hfos;
    organoid.LargeBlockAnalysis.proportion_delta                = proportion_delta;
    organoid.LargeBlockAnalysis.proportion_theta                = proportion_theta;
    organoid.LargeBlockAnalysis.proportion_sgamma               = proportion_sgamma;
    organoid.LargeBlockAnalysis.proportion_fgamma               = proportion_fgamma;
    organoid.LargeBlockAnalysis.proportion_hfos                 = proportion_hfos;
    
    clear stds_largeblocks  coastline_largeblocks norm_coastline_largeblocks ...
        bandpow_all bandpow_delta bandpow_theta bandpow_sgamma bandpow_fgamma bandpow_hfos ...
        proportion_delta  proportion_theta proportion_sgamma proportion_fgamma proportion_hfos

    %% Run analysis of post-KCl perturbations only if the perturbations exist in the recording!
    
    if exist('kclTimeCol','var') && ~isempty(kclTimeCol) && (~strcmpi(raw{org,kclTimeCol},'') && ~isnan(raw{org,kclTimeCol}))
        
        % Normalise std to the std of the total baseline time
        baseline_min = round(raw{org,kclTimeCol}./60); % round it just in case it is expected to be an integer...
        baselinestd = std(data(1:round(baseline_min*60*fs)),1);
        normstds = stds_smallblocks./baselinestd;
        
        % Check if normSTDs go over the thresholds (i.e. >2x, 3x, 4x, 5x... normstd)
		std_thresh = 2:8;
        [std_thresholds] = findNormStdevCrossings(normstds,std_thresh,4);
        
        % Some std seg stuff
        normstd_max = NaN(4,1);
        normstd_avg = NaN(4,1);
        end_of_data_flag = 0; % let's use this to check if we reach the end of the data before expected (i.e. the recording is shorter than expected)
        for segment = 1:4
            seg_ind_end = (segment)*30*60/main_params.block_length_sec;
            if seg_ind_end > length(normstds)
                end_of_data_flag = 1;
                seg_ind_end = length(normstds);
            end
            normstd_max(segment) = max(normstds(1+((segment-1)*30*60/main_params.block_length_sec):seg_ind_end));
            normstd_avg(segment) = nanmean(normstds(1+((segment-1)*30*60/main_params.block_length_sec):seg_ind_end));
            if end_of_data_flag
                break
            end
        end
        clear end_of_data_flag segment

        organoid.SmallBlockAnalysis.perturbation_present = 1;
        organoid.SmallBlockAnalysis.normalised_small_block_analysis_done = 1;
        organoid.SmallBlockAnalysis.stdevs_norm = normstds;
        organoid.SmallBlockAnalysis.std_threshold_inds = std_thresholds;
        organoid.SmallBlockAnalysis.std_threshold_xvals = std_thresh;
        organoid.SmallBlockAnalysis.std_norm_max = normstd_max;
        organoid.SmallBlockAnalysis.std_norm_avg = normstd_avg;
        
    end
    
    

    

    %% SEGMENT Analysis

    length_segment_samples = round(fs*60*segment_length_min);
    if length_segment_samples <= length(data) % this means we have AT LEAST two segments
        
        % construct segments matrix - and then truncate if it exceeds the length of the recording
        % the logic here is really ugly - I think there's a neat way of doing it in a for loop but I can't test this at the moment.
        segments = [1,length_segment_samples;length_segment_samples+1,2*length_segment_samples;2*length_segment_samples+1,3*length_segment_samples;3*length_segment_samples+1,4*length_segment_samples];
        if segments(4,2) > length(data)
            if (segments(4,1)>=length(data)) % uh oh if this isn't true - it means we've lost a whole segment (or the data isn't long enough for the last segment)
                
                segments(4,:) = []; % delete the fourth segment indexer as this doesn't exist
                if(segments(3,2)>length(data))
                    if segments(3,1)>=length(data) % uhoh - there's no third segment either!
                        
                        segments(3,:) = [];
                        if segments(2,2)>length(data)
                            assert(segments(2,1)>length(data)); % this really should be true as we shouldn't be in this top level if-condition unless length(data)>length_segment_samples
                            segments(2,2) = length(data);
                        end
                        
                    else
                        segments(3,2) = length(data);
                    end
                end
                
            else
                segments(4,2) = length(data);
            end
        end
        
        [bandpow_all,~,~] = bandPowerAnalysisSegments(data(:,1),fs,segments);
        [~,bandpow_delta,proportion_delta] = bandPowerAnalysisSegments(data(:,1),fs,segments,main_params.deltaBand,params);
        [~,bandpow_theta,proportion_theta] = bandPowerAnalysisSegments(data(:,1),fs,segments,main_params.thetaBand,params);
        [~,bandpow_sgamma,proportion_sgamma] = bandPowerAnalysisSegments(data(:,1),fs,segments,[25 50],params);
        [~,bandpow_fgamma,proportion_fgamma] = bandPowerAnalysisSegments(data(:,1),fs,segments,[50 100],params);
        [~,bandpow_hfos,proportion_hfos] = bandPowerAnalysisSegments(data(:,1),fs,segments,[100 250],params);

        normpower = bandpow_all./bandpow_all(1);

        organoid.SegmentAnalysis.segment_analysis_done = 1;
        organoid.SegmentAnalysis.segments = segments;
        organoid.SegmentAnalysis.BandpowerAll = bandpow_all;
        organoid.SegmentAnalysis.bandpow_delta = bandpow_delta;
        organoid.SegmentAnalysis.bandpow_theta = bandpow_theta;
        organoid.SegmentAnalysis.bandpow_sgamma = bandpow_sgamma;
        organoid.SegmentAnalysis.bandpow_fgamma = bandpow_fgamma;
        organoid.SegmentAnalysis.bandpow_hfos = bandpow_hfos;
        organoid.SegmentAnalysis.proportion_delta = proportion_delta;
        organoid.SegmentAnalysis.proportion_theta = proportion_theta;
        organoid.SegmentAnalysis.proportion_sgamma = proportion_sgamma;
        organoid.SegmentAnalysis.proportion_fgamma = proportion_fgamma;
        organoid.SegmentAnalysis.proportion_hfos = proportion_hfos;
        organoid.SegmentAnalysis.normPowerAll = normpower;

        % make overall vectors for the orgs instead of saving
        allpower(org,:) = bandpow_all;
        thetapower(org,:) = bandpow_theta;
        sgammapower(org,:) = bandpow_sgamma;
        fgammapower(org,:) = bandpow_fgamma;
    
    else 
        organoid.SegmentAnalysis.segment_analysis_done = 0;
    end
    
    
    %% Entropy analysis
    %{
    
    % I think this can be omitted for now, as shortBlockAnalysis computed
    entropy.  However it does NOT compute normalised entropy, whereas the function below does, and this
    function is probably better in the long run...
    
    [entropy] = organoidEntropy(data,fs);
    organoid.EntropicAnalysis = entropy;
    %}
    
    %% Save the organoid struct and write its location to the xlsfile
    [fpath_excel,~,~] = fileparts(xlsfile);
    save_near_excel_sheet = 1;
    if ~strcmpi(fpath_excel,fpath) && save_near_excel_sheet
        fpath = fpath_excel; 
    end
    orgstructsavefile = sprintf('%s\\%s.mat',fpath,raw{org,orgCol});
    save(orgstructsavefile,'organoid')

    %% Plot organoid summary figure and save

    orgsummarysavefile = sprintf('%s\\Analysis\\%s.png',fpath,raw{org,orgCol});
    
    plotOrgSummaryFigure(organoid,orgsummarysavefile)
    xlswrite(xlsfile,{orgstructsavefile},'Recordings',sprintf('Q%d',org)) % should really convert saveFileCol to Q
    clear fpath fname data organoid orgsummarysavefile

end

[fpath,~,~] = fileparts(xlsfile);
xlsdatafile = sprintf('%s\\OrganoidData.xlsx',fpath);
writeOrganoidStructToExcel(xlsfile,xlsdatafile)

keyboard;