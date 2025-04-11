function [d_ds,fs] = downsampleData(d,si,desired_rate)

%%% FUNCTION downsampleData
%%% si in us - TO DO: make unit specifiable
%%% desired_rate in Hz

if ~exist('si','var')
    error('No sampling interval specified.  Cannot downsample.\n')
end

if ~exist('desired_rate','var')
    fprintf('Warning: No desired rate given!  Using default value of 500Hz\n')
    if si>(1000000/500)
        fprintf('Warning: Data sampling rate <500Hz.  Returning raw trace... \n')
        d_downsampled = d;
        fs = 1/si;
        return;
    else desired_rate = 500;
    end
end

if desired_rate > 1000000/si
    fprintf('Warning: Desired rate greater than actual rate.  Returning raw trace...\n ')
    d_downsampled = d;
    fs = 1/si;
    return;
end

%% Orient and detrend data before anything else

if size(d,1)<size(d,2)
    fprintf('Warning: More datapoints in rows (channels) than columns (datastream)\nI assume there are fewer channels than datapoints - flipping the matrix')
    d = d';
end
d = detrend(d); % detrend operates over COLUMNS of d if d is a matrix

%% Downsample the timeseries to desired rate

%lowPassFilt = designfilt('lowpassfir','PassbandFrequency',500,'StopbandFrequency',1000,'SampleRate',10000);%00./h.si);
%d_filt = filtfilt(lowPassFilt,d);
%{
bandPassFilt = designfilt('bandpassfir','StopbandFrequency1',0.01,'PassbandFrequency1',0.1,'PassbandFrequency2',1000,'StopbandFrequency2',1100,'SampleRate',10000);%00./h.si);
figure; fvtool(bandPassFilt);
d_filt = filtfilt(bandPassFilt,d);
d_downsampled = downsample(d_filt,20);
%}

fprintf('... Downsampling to target %d Hz...\n',desired_rate)

si_desired = 1000000/desired_rate;
f1 = factor(floor(si_desired/ si));
f2 = factor(ceil(si_desired/si));
if length(f1)>length(f2)
    factors = f1;
elseif length(f2)>length(f1)
    factors = f2;
else factors = f1;
end
clear f1 f2

% Alternative way to downsample in MATLAB that is likely quicker and more
% valid:
%d_downsampled = decimate(d,downsample1);
%d_downsampled = decimate(d_downsampled,downsample2);

newsi = repmat(si,size(d,2),1);
for ch = 1:size(d,2)
    d_downsampled = d(:,ch);
    for f = 1:length(factors)
        d_downsampled = decimate(d_downsampled,factors(f));
        newsi(ch) = newsi(ch)*factors(f);
    end
    d_ds(:,ch) = d_downsampled;
    clear d_downsampled
end

si = newsi(1); % These should be the same... test it with eps?
fs = 1000000/si;

fprintf('... Downsampling complete.  Final sampling rate: %d Hz\n', fs);

clear d d_filt
