function [data,params] = preprocessLFP(data,si,params)

%%% FUNCTION preprocessLFP
%%% The aim of this function is to do some of the common early
%%% preprocessing stuff - e.g. detrend, highpass filter 1Hz, bandstop filter 50Hz

if ~exist('params','var')
    params = struct;
end
if ~isfield(params,'detrend')
    params.detrend = 1;
end
if ~isfield(params,'lowpass300Hz')
    params.lowpass300Hz = 1;
end
if ~isfield(params,'highpassCutoff')
    params.lowpassCutoff = 250; % in Hz
end
if ~isfield(params,'highpass1Hz')
    params.highpass1Hz = 1;
end
if ~isfield(params,'highpassCutoff')
    params.highpassCutoff = 1; % in Hz
end
if ~isfield(params,'bandstop50Hz')
    params.bandstop50Hz = 1; % in Hz
end

if ~isfield(params,'save_filter_figures')
    params.save_filter_figures = 0;
elseif ~isfield(params,'save_fpath') || ~exist(params.save_fpath,'dir')
    % else check the filepath exists in params and is real - otherwise don't save figs and warn this is the case
    params.save_filter_figures = 0;
    warning('Organoid:preprocessLFP:noSavePath','No valid save filepath for preprocessing/filters given.  Not saving preprocessing figure...');
end

%% Detrend

if params.detrend
    data = detrend(data); % nice and easy.  Could add in optional n for degree of polynomial fit in future
end


%% Highpass @ 1Hz (or as specified)

if params.highpass1Hz

    fprintf('... Highpass filtering %d Hz\n',params.highpassCutoff)

    highpass1 = designfilt('highpassiir','FilterOrder',8, ...
                 'PassbandFrequency',params.highpassCutoff,'PassbandRipple',0.2, ...
                 'SampleRate',1000000 / si);
    if ~isstable(highpass1)
        keyboard;
    end
    data = filtfilt(highpass1,data);

end
%fvtool(highpass1)

%% Lowpass @ 275Hz

if params.lowpass300Hz

    fprintf('... Lowpass filtering %d Hz\n',params.lowpassCutoff)

    lowpass1 = designfilt('lowpassiir','FilterOrder',8, ...
                 'PassbandFrequency',params.lowpassCutoff,'PassbandRipple',0.2, ...
                 'SampleRate',1000000 / si);
    if ~isstable(lowpass1)
        keyboard;
    end
    data = filtfilt(lowpass1,data);

end

%% Bandstop @ 50Hz + harmonics

if params.bandstop50Hz

    fprintf('... Notch filtering 50Hz + harmonics\n')
    
    notchfilt50 = designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',48,'HalfPowerFrequency2',52,'SampleRate',1000000 / si);
    if ~isstable(notchfilt50)
        keyboard;
    end
    data = filtfilt(notchfilt50,data);
    
    notchfilt100 = designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',98,'HalfPowerFrequency2',102,'SampleRate',1000000 / si);
    if ~isstable(notchfilt100)
        keyboard;
    end
    data = filtfilt(notchfilt100,data);
    
    notchfilt150 = designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',148,'HalfPowerFrequency2',152,'SampleRate',1000000 / si);
    if ~isstable(notchfilt150)
        keyboard;
    end
    data = filtfilt(notchfilt150,data);


end