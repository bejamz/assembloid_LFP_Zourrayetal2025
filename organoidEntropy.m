function [entropy] = organoidEntropy(data,fs,params)

%%% params.
%%%       .segment_length_sec
%%%       .segment_bounds - Do not specify at same time as     segment_length_sec.

if ~exist('params','var') || ~isstruct(params)
    params = struct;
end
if isfield(params,'segment_bounds_sec')
    assert(params.segment_bounds_sec(end)<=(length(data)/60))
    segment_bounds = params.segment_bounds_sec;
    %%% TO DO - check segment bounds is nx2 in size!
elseif isfield(params,'segment_length_sec')
    segment_length_sec = params.segment_length_sec;
    bound_limits = [0:segment_length_sec:(length(data)-1)/fs];
    bound_limits = bound_limits(:);
    segment_bounds = [bound_limits(1:end-1),bound_limits(2:end)-1];
    clear bound_limits
else 
    segment_bounds = [0,(length(data)-1)/fs]; % TimeLimits as an input to pentropy expects bounds within this range [0,length_sec];
end
    
n_segments = size(segment_bounds,1);
entropies = NaN(n_segments,1);
norm_entropies = NaN(n_segments,1);

for seg = 1:n_segments
    
    entropies(seg) = pentropy(data,fs,'TimeLimits',segment_bounds(n_segments,:),'Instantaneous',false,'Scaled',false); % NOTE - TimeLimits [t1 t2] is in seconds if fs is given in Hz
    norm_entropies(seg) = pentropy(data,fs,'TimeLimits',segment_bounds(n_segments,:),'Instantaneous',false,'Scaled',true);
    
end

time_points_sec = nanmean(segment_bounds,2);

entropy = struct;
entropy.time_points_sec = time_points_sec;
entropy.entropies = entropies;
entropy.norm_entropies = norm_entropies;