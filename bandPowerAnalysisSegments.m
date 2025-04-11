function [bandpow_all,bandpow_freqband,bandpow_proportion] = bandPowerAnalysisSegments(data,fs,segments,freqrange,params)

%%% Does: STDEV, Coastline, 
%%% Segments = nx2 matrix of start and end inds of interest
%%% Freqrange should be [a b] where b > a

if ~exist('params','var')
    params = struct;
end
if ~isfield(params,'Method')
    params.Method = 'pwelch';
end


%% check if block_length is in sec and convert to min if so


%% Spectral density calculations for longer blocks

%%% This ignores the concept of baseline/trial segments.  
%%% It will need to be rewritten to allow for these.

fprintf('... Short block analysis\n')

bandpow_all = NaN(size(segments,1),1);
bandpow_freqband = NaN(size(segments,1),1);
bandpow_proportion = NaN(size(segments,1),1);

fprintf('... Progress: 000%%')
for block = 1:size(segments,1)
    
    block_start = segments(block,1);
    block_end   = segments(block,2);

    assert(block_start<block_end);
    assert(block_end<=length(data))
    assert(block_start<=length(data));

    d_block = data(block_start:block_end);

    
    fprintf('\b\b\b\b')
    fprintf('%03d%%',floor(100*block/size(segments,1)))

    bandpow_all(block) = bandpower(d_block,fs,[0 floor(fs/2)]); % get total power up to (below) nyquist freq - for some reason this throws an error if you put in the real nyquist
    if exist('freqrange','var')
        bandpow_freqband(block) = bandpower(d_block,fs,freqrange);
        bandpow_proportion(block) = bandpow_freqband(block)./bandpow_all(block);
    end


    clear d_block block_start block_end 

end
fprintf('\n')




