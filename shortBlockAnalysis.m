function [stds_smallblocks,norm_coastline_smallblocks,coastline_smallblocks,rms_smallblocks,spectralstuff] = shortBlockAnalysis(data,fs,block_length_sec,params)

%%% Does: STDEV, Coastline, Norm Coastline, RMS, Energy
%%% Block length is by default in seconds

if ~exist('params','var')
    params = struct;
end


%% check if block_length is in min and convert to sec if so


%% SD analyses for short blocks

%%% This ignores the concept of baseline/trial segments.  
%%% It will need to be rewritten to allow for these.

fprintf('... Short block analysis\n')

if ~exist('block_length_sec','var')
    block_length_sec = 10;
end

n_samples_per_smallblock = ceil(block_length_sec*fs);
n_smallblocks = ceil((length(data))/n_samples_per_smallblock); %ceil instead of floor for the moment; probably want to reverse this eventually + take blocks from the end backwards.

%% Preallocate vectors

stds_smallblocks = NaN(n_smallblocks,1);
rms_smallblocks = NaN(n_smallblocks,1);
norm_coastline_smallblocks = NaN(n_smallblocks,1);
coastline_smallblocks = NaN(n_smallblocks,1);
medianfreq_smallblocks = NaN(n_smallblocks,1);
meanfreq_smallblocks = NaN(n_smallblocks,1);
energy_smallblocks = NaN(n_smallblocks,1);
entropy_smallblocks = NaN(n_smallblocks,1);
kurtosis_smallblocks = NaN(n_smallblocks,1);
skewness_smallblocks = NaN(n_smallblocks,1);
crest_smallblocks = NaN(n_smallblocks,1);
flatness_smallblocks = NaN(n_smallblocks,1);

% Baseline average stds - nb averaging stds means taking the root of the
% mean of the variances!

%% Begin main loop!

fprintf('... Progress: 000%%')
for block = 1:n_smallblocks
    
    block_start = 1+((block-1).*n_samples_per_smallblock);
    block_end = block_start + n_samples_per_smallblock - 1;
    % Extract this block indices

    if(block_end>length(data))
        block_end = length(data);
    end
    assert(block_start<=length(data));

    d_block = data(block_start:block_end);

    % raw signal stdev
    stds_smallblocks(block) = std(d_block);
    coastline_smallblocks(block) = nansum(abs(diff(d_block))); % in units of voltage of data
    norm_coastline_smallblocks(block) = coastline_smallblocks(block)./(length(d_block)-1); % in units of [V unit/sample]
    
    medianfreq_smallblocks(block) = medfreq(d_block,fs);
    meanfreq_smallblocks(block) = meanfreq(d_block,fs);
    
    % rms and energy
    rms_smallblocks(block) = rms(d_block);
    energy_smallblocks(block) = rms_smallblocks^2 * length(d_block);
    
    % spectral stuff!
    entropy_smallblocks(block) = pentropy(d_block,fs);
    kurtosis_smallblocks(block) = spectralKurtosis(d_block,fs);
    skewness_smallblocks(block) = spectralSkewness(d_block,fs);
    crest_smallblocks(block) = spectralCrest(d_block,fs);
    flatness_smallblocks(block) = spectralFlatness(d_block,fs);
    
    fprintf('\b\b\b\b')
    fprintf('%03d%%',floor(100*block/n_smallblocks))

    clear d_block block_start block_end 

end

spectralstuff = struct;
spectralstuff.signal_energy = energy_smallblocks;
spectralstuff.signal_entropy = entropy_smallblocks;
spectralstuff.spectral_kurtosis = kurtosis_smallblocks;
spectralstuff.spectral_skewness = skewness_smallblocks;
spectralstuff.spectral_crest = crest_smallblocks;
spectralstuff.spectral_flatness = flatness_smallblocks;

fprintf('\n')
