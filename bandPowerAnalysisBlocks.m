function [bandpow_all,bandpow_freqband,bandpow_proportion,block_length_min] = bandPowerAnalysisBlocks(data,fs,block_length_min,freqrange,params)

%%% Does: STDEV, Coastline, 
%%% Block length is by default in seconds
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

if ~exist('block_length_min','var')
    block_length_min = 10;
end

n_samples_per_block = ceil(block_length_min*fs*60);
n_blocks = ceil((length(data))/n_samples_per_block); %ceil instead of floor for the moment; probably want to reverse this eventually + take blocks from the end backwards.


bandpow_all = NaN(n_blocks,1);
bandpow_freqband = NaN(n_blocks,1);
bandpow_proportion = NaN(n_blocks,1);

fprintf('... Progress: 000%%')
for block = 1:n_blocks
    
    block_start = 1+((block-1).*n_samples_per_block);
    block_end = block_start + n_samples_per_block - 1;
    % Extract this block indices

    if(block_end>length(data))
        %fprintf('Note: this block (%d of %d) truncated by %d elements \n',block,n_blocks,block_end-length(data));  % this won't work with the other fprintfs in the loop...
        block_end = length(data);
    end
    assert(block_start<=length(data));

    d_block = data(block_start:block_end);
    n_samples = length(d_block);
    
    fprintf('\b\b\b\b')
    fprintf('%03d%%',floor(100*block/n_blocks))
    
    if logical(mod(n_samples,2)) % we're setting fs_range to be strictly <= max FREQRANGE specified in the bandpower documentation
        fs_range = [0 floor(fs*(n_samples-1)/(2*n_samples))];
    else 
        fs_range = [0 floor(fs/2)];
    end
    
    bandpow_all(block) = bandpower(d_block,fs,fs_range); % get total power up to (below) nyquist freq - for some reason this throws an error if you put in the real nyquist
    
    if exist('freqrange','var')
        bandpow_freqband(block) = bandpower(d_block,fs,freqrange);
        bandpow_proportion(block) = bandpow_freqband(block)./bandpow_all(block);
    end


    clear d_block block_start block_end n_samples

end
fprintf('\n')




