
[num,txt,raw] = xlsread(xlsfile);

%xlsoutputs = 'E:\orgo\OrgoRampingK\OrgoRampingKOutputs.xlsx';

[tmpdir,~,~] = fileparts(xlsfile);
savepath = [tmpdir,'\AnalysisFigures'];

onlyGoodOrgs = 1;
segmentOfInterest_min = [15,20];
ch = 1;

trial_length_min = 120;
col_gp_1 = [1,0,0];
col_gp_2 = [0,0,0];

params = struct;
block_length_sec = 10;
block_length_min = 5;
segment_length_min = 30;

%% TO BE CHANGED!

crispr_inds = zeros(size(raw,1),1);
crispr_inds(53) = 1; crispr_inds(54) = 1;
crispr_inds(51) = 1; crispr_inds(44) = 1;

%% Preallocate arrays

%%% General use arrays
group                       = NaN(size(raw,1),1);
goodOrg                     = NaN(size(raw,1),1);
line                        = cell(size(raw,1),1);

%%% Baseline activity arrays
average_std_baseline        = NaN(size(raw,1),1);
average_power_baseline      = NaN(size(raw,1),1);
delta_power_baseline        = NaN(size(raw,1),1);
theta_power_baseline        = NaN(size(raw,1),1);
beta_power_baseline         = NaN(size(raw,1),1);
sgamma_power_baseline       = NaN(size(raw,1),1);
fgamma_power_baseline       = NaN(size(raw,1),1);
hfo_power_baseline          = NaN(size(raw,1),1);

proportion_delta_baseline   = NaN(size(raw,1),1);
proportion_theta_baseline   = NaN(size(raw,1),1);
proportion_beta_baseline    = NaN(size(raw,1),1);
proportion_sgamma_baseline  = NaN(size(raw,1),1);
proportion_fgamma_baseline  = NaN(size(raw,1),1); 
proportion_hfos_baseline    = NaN(size(raw,1),1);

energy_baseline         = NaN(size(raw,1),1);
skewness_baseline       = NaN(size(raw,1),1);
crest_baseline          = NaN(size(raw,1),1);
entropy_baseline        = NaN(size(raw,1),1);
kurtosis_baseline       = NaN(size(raw,1),1);
flatness_baseline       = NaN(size(raw,1),1);



new_stds = NaN(8,size(raw,1));
allpower = NaN(size(raw,1),4);
thetapower = NaN(size(raw,1),4);
sgammapower = NaN(size(raw,1),4);
fgammapower = NaN(size(raw,1),4);
hfopower = NaN(size(raw,1),4);
normpower = NaN(size(raw,1),4);

proportion_delta        = NaN(size(raw,1),480);
proportion_theta        = NaN(size(raw,1),480);
proportion_beta         = NaN(size(raw,1),480);
proportion_sgamma       = NaN(size(raw,1),480);
proportion_fgamma       = NaN(size(raw,1),480);
proportion_hfos         = NaN(size(raw,1),480);

all_kurtosis            = NaN(size(raw,1),480);
all_skewness            = NaN(size(raw,1),480);
all_energy              = NaN(size(raw,1),480);
all_entropy             = NaN(size(raw,1),480);
all_crest               = NaN(size(raw,1),480);
all_flatness            = NaN(size(raw,1),480);

psd_foldchange_median   = NaN(size(raw,1),1);
psd_foldchange_lowfreq   = NaN(size(raw,1),1);
psd_foldchange_verylowfreq = NaN(size(raw,1),1);
psd_foldchange_15mins = NaN(size(raw,1),1);
psd_all = NaN(size(raw,1),325-8+1);

baseline_std = NaN(size(raw,1),1);
std_segOfInterest = NaN(size(raw,1),1);
std_minsAfterKcl = NaN(size(raw,1),1);

powerspectra = NaN(size(raw,1),50);
entropy = NaN(size(raw,1),1);
norm_entropy = NaN(size(raw,1),1);

criticality_a = NaN(10,2,size(raw,1));
criticality_a_half = NaN(10,2,size(raw,1));
criticality_a_quarter = NaN(10,2,size(raw,1));
criticality_a_lagone = NaN(size(raw,1),480);
criticality_a_norm_baseline = NaN(size(raw,1),1);
criticality_a_norm_4std = NaN(size(raw,1),1);
criticality_a_norm_1std = NaN(size(raw,1),1);

normstds_median = NaN(480,size(raw,1));



%% Get meta info about the analyses - notable small block length (sec)

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
            lineCol = col;
        case 'date'
            dateCol = col;
        case 'batch' 
            batchCol = col;
        case 'age'
            ageCol = col;
        case 'goodorgmanual' % use 'gabri' for old classification gabri did (does not include batch 3 at the moment)
            goodOrgCol = col;
    end
end

tmp_org = 3; %This will be the first organoid in the excel list
load(sprintf('%s\\%s.mat',raw{tmp_org,fpathCol},raw{tmp_org,orgCol}));
small_block_sec = organoid.SmallBlockAnalysis.Meta.SmallBlockLengthSec;
n_blocks = trial_length_min*60/small_block_sec;
normstds_all = NaN(n_blocks,size(raw,1));
normcoastline_all = NaN(n_blocks,size(raw,1));

energy = NaN(size(raw,1),n_blocks);

clear organoid tmp_org fname 

normstds_all = NaN(n_blocks,size(raw,1));
normcoastline_all = NaN(n_blocks,size(raw,1));
std_thresh_all = NaN(size(raw,1),7);

all_power_15min = NaN(size(raw,1),8);
all_std_15min = NaN(size(raw,1),8);

fprintf('%d organoids found...\n',size(raw,1)-2);

for org = 3:size(raw,1)
    
    fprintf('... running org no. %d ...',org-2);
    
    if isnan(raw{org,orgCol})
        fprintf('skipping...\n')
        continue;
    end
    
    if ~raw{org,includeCol} % if not due to be included in the dataset
        fprintf('skipping...\n')
        continue;
    end
    
    group(org) = raw{org,groupCol};
    goodOrg(org) = raw{org,goodOrgCol};
    line{org}    = raw{org,lineCol};
    
    assert(exist(raw{org,saveFileCol},'file')==2)
    load(raw{org,saveFileCol})
    
    % Plot the actual trace
    %{
    clf;
    if raw{org,7} == 2
        plot([1:size(organoid.Data.data,1)]/(organoid.Data.Fs*60),organoid.Data.data(:,1),'Color','k')
        xlim([0 120]); xlabel('Time /min')
        ylim([nanmean(organoid.Data.data(:,1))-20*nanstd(organoid.Data.data(:,1)), nanmean(organoid.Data.data(:,1))+20*nanstd(organoid.Data.data(:,1))])
        ylabel('Voltage /mV')
        title(sprintf('%s',raw{org,1}))
        saveas(gcf,sprintf('E:\\orgo\\OrgoRampingK\\Figures\\%s.png',raw{org,1}))
    elseif raw{org,7} == 1
        plot([1:size(organoid.Data.data,1)]/(organoid.Data.Fs*60),organoid.Data.data(:,1),'Color','r')
        xlim([0 120]); xlabel('Time /min')
        ylim([nanmean(organoid.Data.data(:,1))-20*nanstd(organoid.Data.data(:,1)), nanmean(organoid.Data.data(:,1))+20*nanstd(organoid.Data.data(:,1))])
        ylabel('Voltage /mV')
        title(sprintf('%s',raw{org,1}))
        saveas(gcf,sprintf('E:\\orgo\\OrgoRampingK\\Figures\\%s.png',raw{org,1}))
    end
    %}
    
    fs = organoid.Data.Fs;
    kcl_time_sec = raw{org,kclTimeCol};
    short_block_length_sec = organoid.SmallBlockAnalysis.Meta.SmallBlockLengthSec;
    
    %% TEMPORARY!!!  if .EDR then fix the scaling issue...
    if strcmpi(organoid.Data.meta.fileType,'EDR')
        disp(' Alert: fixing the EDR scaling!')
        rescaling_factor = (organoid.Data.meta.AD./(organoid.Data.meta.ADCMAX*organoid.Data.meta.YCF2));
        fprintf('... Rescaling Y-axis by %3.2f ... ',1/rescaling_factor);
        organoid.Data.data = organoid.Data.data .* rescaling_factor;
    end
    
    
    if n_blocks > length(organoid.SmallBlockAnalysis.stdevs_norm)
        end_ind = length(organoid.SmallBlockAnalysis.stdevs_norm);
    else 
        end_ind = n_blocks;
    end
    
    %% Extract baseline info
    
    tmp_data = organoid.Data.data(1:round(fs*kcl_time_sec),1);
    average_std_baseline(org)       = nanstd(tmp_data);
    average_power_baseline(org)     = bandpower(tmp_data);
    
    delta_power_baseline(org)    = nanmean(organoid.SmallBlockAnalysis.bandpow_delta(1:floor(kcl_time_sec/short_block_length_sec),1));
    theta_power_baseline(org)    = nanmean(organoid.SmallBlockAnalysis.bandpow_theta(1:floor(kcl_time_sec/short_block_length_sec),1));
    beta_power_baseline(org)     = nanmean(organoid.SmallBlockAnalysis.bandpow_beta(1:floor(kcl_time_sec/short_block_length_sec),1));
    sgamma_power_baseline(org)   = nanmean(organoid.SmallBlockAnalysis.bandpow_sgamma(1:floor(kcl_time_sec/short_block_length_sec),1));
    fgamma_power_baseline(org)   = nanmean(organoid.SmallBlockAnalysis.bandpow_fgamma(1:floor(kcl_time_sec/short_block_length_sec),1));
    hfo_power_baseline(org)      = nanmean(organoid.SmallBlockAnalysis.bandpow_hfos(1:floor(kcl_time_sec/short_block_length_sec),1));

    proportion_delta_baseline(org)   = nanmean(organoid.SmallBlockAnalysis.proportion_delta(1:floor(kcl_time_sec/short_block_length_sec),1));
    proportion_theta_baseline(org)   = nanmean(organoid.SmallBlockAnalysis.proportion_theta(1:floor(kcl_time_sec/short_block_length_sec),1));
    proportion_beta_baseline(org)    = nanmean(organoid.SmallBlockAnalysis.proportion_beta(1:floor(kcl_time_sec/short_block_length_sec),1));
    proportion_sgamma_baseline(org)  = nanmean(organoid.SmallBlockAnalysis.proportion_sgamma(1:floor(kcl_time_sec/short_block_length_sec),1));
    proportion_fgamma_baseline(org)  = nanmean(organoid.SmallBlockAnalysis.proportion_fgamma(1:floor(kcl_time_sec/short_block_length_sec),1));
    proportion_hfos_baseline(org)    = nanmean(organoid.SmallBlockAnalysis.proportion_hfos(1:floor(kcl_time_sec/short_block_length_sec),1));
    
    %{
    %%% NEEDS TO BE COMPUTED LOCALLY UNTIL THE .EDR LOADIN FUNCTION IS FIXED
    energy_baseline(org)         = nanmean(organoid.SmallBlockAnalysis.SpectralMeasurements.signal_energy(1:floor(kcl_time_sec/short_block_length_sec),1));
    skewness_baseline(org)       = nanmean(organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_skewness(1:floor(kcl_time_sec/short_block_length_sec),1));
    crest_baseline(org)          = nanmean(organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_crest(1:floor(kcl_time_sec/short_block_length_sec),1));
    entropy_baseline(org)        = nanmean(organoid.SmallBlockAnalysis.SpectralMeasurements.signal_entropy(1:floor(kcl_time_sec/short_block_length_sec),1));
    kurtosis_baseline(org)       = nanmean(organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_kurtosis(1:floor(kcl_time_sec/short_block_length_sec),1));
    flatness_baseline(org)       = nanmean(organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_flatness(1:floor(kcl_time_sec/short_block_length_sec),1));
    %}
    
    energy_baseline(org)         = rms(tmp_data)^2 * length(tmp_data);
    skewness_baseline(org)       = spectralSkewness(tmp_data,fs,'Window',ones(size(tmp_data)));
    crest_baseline(org)          = spectralCrest(tmp_data,fs,'Window',ones(size(tmp_data)));
    entropy_baseline(org)        = spectralEntropy(tmp_data,fs,'Window',ones(size(tmp_data)));
    kurtosis_baseline(org)       = spectralKurtosis(tmp_data,fs,'Window',ones(size(tmp_data)));
    flatness_baseline(org)       = spectralFlatness(tmp_data,fs,'Window',ones(size(tmp_data)));
    
    clear tmp_data
    
    %% Overall power changes post-KCl (i.e. averaged through time)
    
    %%% 15 minute chunks - STD and Power
    for seg = 1:8 % NB - this is hard coded!  assumes kcl time is 15 mins and recordings are 2 hours
        start_ind = round((fs*60*15*(seg-1))+1);
        end_ind_tmp = round((fs*60*15*(seg)));
        if start_ind>length(organoid.Data.data) || end_ind_tmp > length(organoid.Data.data)
            continue
        end
        data_chunk = organoid.Data.data(start_ind:end_ind_tmp);
        all_power_15min(org,seg) = bandpower(data_chunk);
        all_std_15min(org,seg) = nanstd(data_chunk);
        clear data_chunk start_ind end_ind_tmp seg
    end
    
    %%% Extract fold change in PSD
    %%%% NB: WHY DIDN'T I SAVE THE PSD??  Have to compute it now...
    tmp_data = organoid.Data.data(1:round(fs*kcl_time_sec),1);
    [pxx,fxx] = pwelch(tmp_data,round(fs*5),[],[],fs);
    [pxx2,fxx2] = pwelch(organoid.Data.data(round(fs*kcl_time_sec)+1:end,1),round(fs*5),[],[],fs); % TODO: change "end" to exactly 2 hours (or end if 2 hours is longer)
    [pxx3,fxx3] = pwelch(organoid.Data.data(round(fs*kcl_time_sec)+1:round(fs*kcl_time_sec)+round(fs*900),1),round(fs*5),[],[],fs); % takes 15 minutes after KCl addition
    
    psd_foldchange_median(org) = median(pxx2(8:end-500)./pxx(8:end-500)); % TO DO - change the end-500 to be specific to up to 250Hz
    %%% TO DO - add psd fold change for different bands?
    psd_foldchange_lowfreq(org) = median(pxx2(8:325)./pxx(8:325)); % TODO - 325 is approx 50Hz
    psd_foldchange_verylowfreq(org) = median(pxx2(8:100)./pxx(8:100)); %TODO 100 is approx 15Hz
    
    psd_foldchange_15mins(org) = median(pxx3(8:325)./pxx(8:325));
    
    psd_all(org,:) = pxx2(8:325)./pxx(8:325);
    
    %% Other stuff
    
    
    baseline_min = 15;
    baselinestd = nanstd(organoid.Data.data(1:round(baseline_min*60*fs)),1);
    %normstds = organoid.SmallBlockAnalysis.stdevs_smallblocks./baselinestd;
    %normstds_all(:,org) = organoid.SmallBlockAnalysis.stdevs_smallblocks(1:n_blocks)./baselinestd;
    normstds_all(1:end_ind,org) = organoid.SmallBlockAnalysis.stdevs_norm(1:end_ind);
    
    %for i = 1:8
    %    new_stds(i,org) = std(organoid.Data.data(round(1+(i-1)*15*60*fs):round(i*15*60*fs),ch),1)./baselinestd;
    %end
    
    std_segOfInterest(org) = nanstd(organoid.Data.data([round(fs*60*segmentOfInterest_min(1)):round(fs*60*segmentOfInterest_min(2))],ch),1)./baselinestd;
    baseline_std(org) = baselinestd;
    
    normcoastline = nanmean(organoid.SmallBlockAnalysis.norm_coastline_smallblocks(1:60));
    normnormcoastline = organoid.SmallBlockAnalysis.norm_coastline_smallblocks./normcoastline;
    normcoastline_all(1:end_ind,org) = normnormcoastline(1:end_ind);
     
    energy(org,1:end_ind) = cumsum(organoid.SmallBlockAnalysis.SpectralMeasurements.signal_energy(1:end_ind));
    
    %% Std relative to the KCl in time - if this varies between organoids
    
    kcl_time_sec = raw{org,kclTimeCol};
    std_baseline_kcl = nanstd(organoid.Data.data(1:round(kcl_time_sec*fs)),1);
    std_segs_mins = 5;
    end_ind_kcl = round(kcl_time_sec*fs + std_segs_mins*60*fs);
    std_norm_kcl = nanstd(organoid.Data.data(round(kcl_time_sec*fs)+1:end_ind_kcl),1)/std_baseline_kcl; %Take the 5 mins after KCl addition
    std_minsAfterKcl(org) = std_norm_kcl;
    
    normstds_median(1:end_ind,org) = organoid.SmallBlockAnalysis.stdevs_smallblocks(1:end_ind)./nanmedian(organoid.SmallBlockAnalysis.stdevs_smallblocks(1:round(kcl_time_sec/organoid.SmallBlockAnalysis.Meta.SmallBlockLengthSec)));
    
    
    %% Extract bandpowers
    thetapower(org,1:4) = organoid.SegmentAnalysis.bandpow_theta;
    sgammapower(org,1:4) = organoid.SegmentAnalysis.bandpow_sgamma;
    fgammapower(org,1:4) = organoid.SegmentAnalysis.bandpow_fgamma;
    hfopower(org,1:4) = organoid.SegmentAnalysis.bandpow_hfos;
    allpower(org,1:4) = organoid.SegmentAnalysis.BandpowerAll;
    normpower(org,1:4) = organoid.SegmentAnalysis.normPowerAll;
    
    %if raw{org,7} == 1
    %    plot(1:length(normstds),normstds,'Color',[0.5,0.5,0.5]);
    %elseif raw{org,7} ==2
    %    plot(1:length(normstds),normstds,'Color',[1,0,0]);
    %end
    
    
    tmp_ind = min([size(proportion_delta,2),length(organoid.SmallBlockAnalysis.proportion_delta)]);
    proportion_delta(org,1:tmp_ind) = organoid.SmallBlockAnalysis.proportion_delta(1:tmp_ind);
    tmp_ind = min([size(proportion_theta,2),length(organoid.SmallBlockAnalysis.proportion_theta)]);
    proportion_theta(org,1:tmp_ind) = organoid.SmallBlockAnalysis.proportion_theta(1:tmp_ind);
    tmp_ind = min([size(proportion_beta,2),length(organoid.SmallBlockAnalysis.proportion_beta)]);
    proportion_beta(org,1:tmp_ind) = organoid.SmallBlockAnalysis.proportion_beta(1:tmp_ind);
    tmp_ind = min([size(proportion_sgamma,2),length(organoid.SmallBlockAnalysis.proportion_sgamma)]);
    proportion_sgamma(org,1:tmp_ind) = organoid.SmallBlockAnalysis.proportion_sgamma(1:tmp_ind);
    tmp_ind = min([size(proportion_fgamma,2),length(organoid.SmallBlockAnalysis.proportion_fgamma)]);
    proportion_fgamma(org,1:tmp_ind) = organoid.SmallBlockAnalysis.proportion_fgamma(1:tmp_ind);
    tmp_ind = min([size(proportion_hfos,2),length(organoid.SmallBlockAnalysis.proportion_hfos)]);
    proportion_hfos(org,1:tmp_ind) = organoid.SmallBlockAnalysis.proportion_hfos(1:tmp_ind);
    
    tmp_ind = min([size(all_kurtosis,2),length(organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_kurtosis)]);
    all_kurtosis(org,1:tmp_ind) = organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_kurtosis(1:tmp_ind);
    tmp_ind = min([size(all_skewness,2),length(organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_skewness)]);
    all_skewness(org,1:tmp_ind) = organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_skewness(1:tmp_ind);
    tmp_ind = min([size(all_energy,2),length(organoid.SmallBlockAnalysis.SpectralMeasurements.signal_energy)]);
    all_energy(org,1:tmp_ind) = organoid.SmallBlockAnalysis.SpectralMeasurements.signal_energy(1:tmp_ind);
    tmp_ind = min([size(all_entropy,2),length(organoid.SmallBlockAnalysis.SpectralMeasurements.signal_entropy)]);
    all_entropy(org,1:tmp_ind) = organoid.SmallBlockAnalysis.SpectralMeasurements.signal_entropy(1:tmp_ind);
    tmp_ind = min([size(all_crest,2),length(organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_crest)]);
    all_crest(org,1:tmp_ind) = organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_crest(1:tmp_ind);
    tmp_ind = min([size(all_flatness,2),length(organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_flatness)]);
    all_flatness(org,1:tmp_ind) = organoid.SmallBlockAnalysis.SpectralMeasurements.spectral_flatness(1:tmp_ind);
    
    clear tmp_ind
    
    std_thresh_all(org,:) = organoid.SmallBlockAnalysis.std_threshold_inds(:);

    %% Entropy and criticality stuff
    
    %%% TO DO - will need to change this when I fix the entropy analysis to
    %%% actually work instead of taking the whole signal!
    if isfield(organoid,'EntropyAnalysis')
        norm_entropy(org)   = organoid.EntropyAnalysis.norm_entropies;
        entropy(org)        = organoid.EntropyAnalysis.entropies;
    end
    
    if isfield(organoid,'Criticality')
        criticality_a_half(:,:,org) = organoid.Criticality.AutocorrWidth.acfw_values_half;
        criticality_a_quarter(:,:,org) = organoid.Criticality.AutocorrWidth.acfw_values_quarter;
        criticality_a(:,:,org) = organoid.Criticality.AutocorrWidth.a_values;
        end_ind = min([480,length(organoid.Criticality.AutocorrWidth.acfw_half)]);
        tmp = organoid.Criticality.AutocorrWidth.acfw_half(1:end_ind)./organoid.Criticality.AutocorrWidth.median_acfw_baseline;
        criticality_a_norm_baseline(org) = nanmean(tmp(1:60));
        inds = find(organoid.SmallBlockAnalysis.stdevs_norm(1:end_ind)>4);
        inds(inds<61) = [];
        criticality_a_norm_4std(org) = nanmean(tmp(inds));
        inds = find(organoid.SmallBlockAnalysis.stdevs_norm(1:end_ind)>0.5 & organoid.SmallBlockAnalysis.stdevs_norm(1:end_ind)<1.5);
        inds(inds<61) = [];
        criticality_a_norm_1std(org) = nanmean(tmp(inds));

        criticality_a_lagone(org,1:end_ind) = organoid.Criticality.AutocorrWidth.acfw_lagone(1:end_ind);
    end
    
    clear end_ind inds
    
    %% keyboard if necessary
    %if org == 19 || org==20
    %    keyboard
    %end

    %% Extract power spectra info
    % to be put in main analysis...

    fprintf('\n')
    
    clear fpath fname data organoid normnormcoastline normcoastline normstds 

end

close all % get rid of any unwanted figures that e.g. pentropy may have dragged up

%% We can change the group variable to set what ones to include
lines_to_include = {'P1','X1','C5'}; % Exclude P3 because of the mutation...
inds_to_include = zeros(size(group));
for i = 1:length(lines_to_include)
    inds_to_include = inds_to_include | strcmpi(line,lines_to_include{i});
end
clear i
group(~inds_to_include) = 0;

%% Clean up the vars with NaNs in them

crispr_inds = crispr_inds(group==1|group==2);
line = line(group==1|group==2);

%%% Clean up baseline vars
average_std_baseline    = average_std_baseline(group==1|group==2,:);
average_power_baseline  = average_power_baseline(group==1|group==2,:);
delta_power_baseline    = delta_power_baseline(group==1|group==2,:);
theta_power_baseline    = theta_power_baseline(group==1|group==2,:);
beta_power_baseline     = beta_power_baseline(group==1|group==2,:);
sgamma_power_baseline   = sgamma_power_baseline(group==1|group==2,:);
fgamma_power_baseline   = fgamma_power_baseline(group==1|group==2,:);
hfo_power_baseline      = hfo_power_baseline(group==1|group==2,:);
proportion_delta_baseline   = proportion_delta_baseline(group==1|group==2,:);
proportion_theta_baseline   = proportion_theta_baseline(group==1|group==2,:);
proportion_beta_baseline    = proportion_beta_baseline(group==1|group==2,:);
proportion_sgamma_baseline  = proportion_sgamma_baseline(group==1|group==2,:);
proportion_fgamma_baseline  = proportion_fgamma_baseline(group==1|group==2,:);
proportion_hfos_baseline    = proportion_hfos_baseline(group==1|group==2,:);
energy_baseline         = energy_baseline(group==1|group==2,:);
skewness_baseline       = skewness_baseline(group==1|group==2,:);
crest_baseline          = crest_baseline(group==1|group==2,:);
entropy_baseline        = entropy_baseline(group==1|group==2,:);
kurtosis_baseline       = kurtosis_baseline(group==1|group==2,:);
flatness_baseline       = flatness_baseline(group==1|group==2,:);

%%% Clean up general power changes
all_power_15min             = all_power_15min(group==1|group==2,:);
all_std_15min               = all_std_15min(group==1|group==2,:);
psd_foldchange_median       = psd_foldchange_median(group==1|group==2);
psd_foldchange_lowfreq      = psd_foldchange_lowfreq(group==1|group==2);
psd_foldchange_verylowfreq  = psd_foldchange_verylowfreq(group==1|group==2);
psd_foldchange_15mins       = psd_foldchange_15mins(group==1|group==2);
psd_all                     = psd_all(group==1|group==2,:);

allpower = allpower(group==1|group==2,:);
normpower = normpower(group==1|group==2,:);
normstds_all = normstds_all(:,group==1|group==2);
normstds_median = normstds_median(:,group==1|group==2);
normcoastline_all = normcoastline_all(:,group==1|group==2);
thetapower = thetapower(group==1|group==2,:);
sgammapower = sgammapower(group==1|group==2,:);
fgammapower = fgammapower(group==1|group==2,:);
hfopower = hfopower(group==1|group==2,:);
%new_stds = new_stds(:,group==1|group==2);
baseline_std = baseline_std(group==1|group==2);
std_thresh_all = std_thresh_all(group==1|group==2,:);
std_minsAfterKcl = (group==1|group==2);

proportion_delta = proportion_delta(group==1|group==2,:);
proportion_theta = proportion_theta(group==1|group==2,:);
proportion_sgamma = proportion_sgamma(group==1|group==2,:);
proportion_fgamma = proportion_fgamma(group==1|group==2,:);
proportion_hfos = proportion_hfos(group==1|group==2,:);

all_kurtosis            = all_kurtosis(group==1|group==2,:);
all_skewness            = all_skewness(group==1|group==2,:);
all_energy              = all_energy(group==1|group==2,:);
all_entropy             = all_entropy(group==1|group==2,:);
all_flatness            = all_flatness(group==1|group==2,:);
all_crest               = all_crest(group==1|group==2,:);

criticality_a_quarter = criticality_a_quarter(:,:,group==1|group==2);
criticality_a_half = criticality_a_half(:,:,group==1|group==2);
criticality_a = criticality_a(:,:,group==1|group==2);
criticality_a_lagone = criticality_a_lagone(group==1|group==2,:);

criticality_a_norm_baseline = criticality_a_norm_baseline(group==1|group==2);
criticality_a_norm_4std = criticality_a_norm_4std(group==1|group==2);
criticality_a_norm_1std = criticality_a_norm_1std(group==1|group==2);

energy = energy(group==1|group==2,:);

goodOrg = goodOrg(group==1|group==2,:);
group = group(~isnan(group));


%% Some of baseline stuff should be done on all organoids, not only the "good" organoids

h1 = figure;

% power stuff

%%% RAW STDEV
cla
ax = subplot(1,2,1);
hold on;
raincloudPlotForPaper(ax,average_std_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,average_std_baseline(group==2),1,0.5,[0,0,0]);
ylim([0 0.5])
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
xlabel('Group')
ylabel('Raw StdDev')
title('Baseline StDev (ylim-[0,0.5])')
ax = subplot(1,2,2);
hold on;
raincloudPlotForPaper(ax,average_std_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,average_std_baseline(group==2),1,0.5,[0,0,0]);
ylim([0 0.05])
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
xlabel('Group')
ylabel('Raw StdDev')
title('Baseline StDev (ylim-[0,0.05])')
set(h1,'Units','Normalized');
set(h1,'Position',[0.3536 0.3454 0.2339 0.5593]);
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineRawStDev,AllOrgs.png',savepath));
[p,h,stats] = ranksum(average_std_baseline(group==1),average_std_baseline(group==2));
clear p h stats

%%% RAW POWER 
clf
ax = gca;
hold on;
raincloudPlotForPaper(ax,log10(average_power_baseline(group==1)),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,log10(average_power_baseline(group==2)),1,0.5,[0,0,0]);
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
ylim([-5 0])
xlabel('Group')
ylabel('Log10(Raw Power)')
title('Baseline power, yscale-log10')
set(h1,'Units','Normalized');
set(h1,'Position',[0.3536 0.3454 0.2339 0.5593]);
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineRawPower,LogScale,AllOrgs.png',savepath));

[p,h,stats] = ranksum(average_power_baseline(group==1),average_power_baseline(group==2))
clear p h stats

%%% PROPORTIONS OF POWER IN FREQUENCY BANDS
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,proportion_delta_baseline(group==1),0,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_delta_baseline(group==2),0.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_theta_baseline(group==1),1,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_theta_baseline(group==2),1.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_beta_baseline(group==1),2,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_beta_baseline(group==2),2.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_sgamma_baseline(group==1),3,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_sgamma_baseline(group==2),3.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_fgamma_baseline(group==1),4,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_fgamma_baseline(group==2),4.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_hfos_baseline(group==1),5,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_hfos_baseline(group==2),5.4,0.2,[0,0,0]);
ylim([0 0.5])
ax.XTick = [0.2,1.2,2.2,3.2,4.2,5.2];
ax.XTickLabel = {'Delta','Theta','Beta','SGamma','FGamma','HFO'};
ylabel('Proportion of power')
xlabel('Frequency band')
set(h1,'Position',[0.1880 0.3769 0.4146 0.4556])
title(sprintf('Mean proportion of power per frequency bands - all orgs'))
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineProportionalPower,AllOrgs.png',savepath));

%%% SPECTRAL ENERGY
clf
ax = gca;
hold on;
raincloudPlotForPaper(ax,log10(energy_baseline(group==1)),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,log10(energy_baseline(group==2)),1,0.5,[0,0,0]);
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
ylim([0 5])
xlabel('Group')
ylabel('Log10(Spectral energy)')
title('Baseline energy, yscale-log10')
set(h1,'Units','Normalized');
set(h1,'Position',[0.3536 0.3454 0.1938 0.5490]);
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineSpectralEnergy,AllOrgs.png',savepath));

%%% SPECTRAL KURTOSIS
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,kurtosis_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,kurtosis_baseline(group==2),1,0.5,[0,0,0]);
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
ylim([1.5 8]) % NOTE!!!! THERE IS ONE ORGANOID WITH KURTOSIS ~90!
xlabel('Group')
ylabel('Spectral kurtosis')
title('Baseline kurtosis, yscale-linear')
set(h1,'Units','Normalized');
set(h1,'Position',[0.3536 0.3454 0.1938 0.5490]);
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineSpectralKurtosis,AllOrgs.png',savepath));

%%% SPECTRAL SKEWNESS
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,skewness_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,skewness_baseline(group==2),1,0.5,[0,0,0]);
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
ylim([-0.5 2.5])
xlabel('Group')
ylabel('Spectral skewness')
title('Baseline skewness, yscale-linear')
set(h1,'Units','Normalized');
set(h1,'Position',[0.3536 0.3454 0.1938 0.5490]);
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineSpectralSkewness,AllOrgs.png',savepath));

%%% SPECTRAL FLATNESS
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,flatness_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,flatness_baseline(group==2),1,0.5,[0,0,0]);
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
ylim([0 0.25])
xlabel('Group')
ylabel('Spectral flatness')
title('Baseline flatness, yscale-linear')
set(h1,'Units','Normalized');
set(h1,'Position',[0.3536 0.3454 0.1938 0.5490]);
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineSpectralFlatness,AllOrgs.png',savepath));

%%% SPECTRAL CREST
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,crest_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,crest_baseline(group==2),1,0.5,[0,0,0]);
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
ylim([0 8000])
xlabel('Group')
ylabel('Spectral crest')
title('Baseline crest, yscale-linear')
set(h1,'Units','Normalized');
set(h1,'Position',[0.3536 0.3454 0.1938 0.5490]);
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineSpectralCrest,AllOrgs.png',savepath));

%%% SPECTRAL ENTROPY
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,entropy_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,entropy_baseline(group==2),1,0.5,[0,0,0]);
ax.XTick = [0,1];
ax.XTickLabel = {'Dravet','Control'};
ylim([0.6 1])
xlabel('Group')
ylabel('Spectral entropy')
title('Baseline entropy, yscale-linear')
set(h1,'Units','Normalized');
set(h1,'Position',[0.3536 0.2472 0.1938 0.6472]);
saveas(h1,sprintf('%s\\Baseline\\AllOrganoids\\BaselineSpectralEntropy,AllOrgs.png',savepath));

close(h1)
%%% TO DO!

%% If interested, truncate to "good" organoids

if onlyGoodOrgs
    goodOrg(isnan(goodOrg))=0;
    goodOrg = logical(goodOrg);
    line = line(goodOrg);
    crispr_inds = crispr_inds(goodOrg);
    
    %%% Baseline vars
    average_std_baseline    = average_std_baseline(goodOrg,:);
    average_power_baseline  = average_power_baseline(goodOrg,:);
    delta_power_baseline    = delta_power_baseline(goodOrg,:);
    theta_power_baseline    = theta_power_baseline(goodOrg,:);
    beta_power_baseline     = beta_power_baseline(goodOrg,:);
    sgamma_power_baseline   = sgamma_power_baseline(goodOrg,:);
    fgamma_power_baseline   = fgamma_power_baseline(goodOrg,:);
    hfo_power_baseline      = hfo_power_baseline(goodOrg,:);
    proportion_delta_baseline   = proportion_delta_baseline(goodOrg,:);
    proportion_theta_baseline   = proportion_theta_baseline(goodOrg,:);
    proportion_beta_baseline    = proportion_beta_baseline(goodOrg,:);
    proportion_sgamma_baseline  = proportion_sgamma_baseline(goodOrg,:);
    proportion_fgamma_baseline  = proportion_fgamma_baseline(goodOrg,:);
    proportion_hfos_baseline    = proportion_hfos_baseline(goodOrg,:);
    energy_baseline         = energy_baseline(goodOrg,:);
    skewness_baseline       = skewness_baseline(goodOrg,:);
    crest_baseline          = crest_baseline(goodOrg,:);
    entropy_baseline        = entropy_baseline(goodOrg,:);
    kurtosis_baseline       = kurtosis_baseline(goodOrg,:);
    flatness_baseline       = flatness_baseline(goodOrg,:);
    
    %%% General changes in power 
    all_power_15min = all_power_15min(goodOrg,:);
    all_std_15min = all_std_15min(goodOrg,:);
    psd_foldchange_median = psd_foldchange_median(goodOrg);
    psd_foldchange_lowfreq = psd_foldchange_lowfreq(goodOrg);
    psd_foldchange_verylowfreq = psd_foldchange_verylowfreq(goodOrg);
    psd_foldchange_15mins = psd_foldchange_15mins(goodOrg);
    psd_all                     = psd_all(goodOrg,:);
    
    allpower = allpower(goodOrg,:);
    normpower = normpower(goodOrg,:);
    normstds_all = normstds_all(:,goodOrg);
    normcoastline_all = normcoastline_all(:,goodOrg);
    normstds_median = normstds_median(:,goodOrg);
    thetapower = thetapower(goodOrg,:);
    sgammapower = sgammapower(goodOrg,:);
    fgammapower = fgammapower(goodOrg,:);
    hfopower = hfopower(goodOrg,:);
    %new_stds = new_stds(:,goodOrg);
    baseline_std = baseline_std(goodOrg);
    std_thresh_all = std_thresh_all(goodOrg,:);
    std_segOfInterest = std_segOfInterest(goodOrg);
    group = group(goodOrg);
    std_minsAfterKcl = std_minsAfterKcl(goodOrg);
    
    proportion_delta = proportion_delta(goodOrg,:);
    proportion_theta = proportion_theta(goodOrg,:);
    proportion_sgamma = proportion_sgamma(goodOrg,:);
    proportion_fgamma = proportion_fgamma(goodOrg,:);
    proportion_hfos = proportion_hfos(goodOrg,:);
    
    
    all_kurtosis            = all_kurtosis(goodOrg,:);
    all_skewness            = all_skewness(goodOrg,:);
    all_energy              = all_energy(goodOrg,:);
    all_entropy             = all_entropy(goodOrg,:);
    all_flatness            = all_flatness(goodOrg,:);
    all_crest               = all_crest(goodOrg,:);
    
    criticality_a_quarter = criticality_a_quarter(:,:,goodOrg);
    criticality_a_half = criticality_a_half(:,:,goodOrg);
    criticality_a = criticality_a(:,:,goodOrg);
    criticality_a_lagone = criticality_a_lagone(goodOrg,:);
    
    criticality_a_norm_baseline = criticality_a_norm_baseline(goodOrg);
    criticality_a_norm_4std = criticality_a_norm_4std(goodOrg);
    criticality_a_norm_1std = criticality_a_norm_1std(goodOrg);
    
    energy = energy(goodOrg,:);
    
    goodOrg = goodOrg(goodOrg);
end


%% Make boolean indices for each line
% Possible lines = P1, P2, P3, C5, X1

num_lines = length(unique(line));
assert(num_lines <= 5); % the below code needs rewriting if we add new lines

p1_inds = logical(zeros(length(line),1));
p2_inds = logical(zeros(length(line),1));
p3_inds = logical(zeros(length(line),1));
c5_inds = logical(zeros(length(line),1));
x1_inds = logical(zeros(length(line),1));
for tmp_ind = 1:length(line)
    switch line{tmp_ind}
        case 'P1'
            p1_inds(tmp_ind) = 1;
        case 'P2'
            p2_inds(tmp_ind) = 1;
        case 'P3'
            p3_inds(tmp_ind) = 1;
        case 'C5'
            c5_inds(tmp_ind) = 1;
        case 'X1'
            x1_inds(tmp_ind) = 1;
    end
end
clear tmp_ind 

%% Baseline activity
% up to the kcl time

h1 = figure;

% power stuff
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,average_std_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,average_std_baseline(group==2),1,0.5,[0,0,0]);
[p,h,stats] = ranksum(average_std_baseline(group==1),average_std_baseline(group==2))
clear p h stats
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,average_power_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,average_power_baseline(group==2),1,0.5,[0,0,0]);
[p,h,stats] = ranksum(average_power_baseline(group==1),average_power_baseline(group==2))
clear p h stats

% Proportional plots
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,proportion_delta_baseline(group==1),0,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_delta_baseline(group==2),0.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_theta_baseline(group==1),1,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_theta_baseline(group==2),1.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_beta_baseline(group==1),2,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_beta_baseline(group==2),2.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_sgamma_baseline(group==1),3,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_sgamma_baseline(group==2),3.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_fgamma_baseline(group==1),4,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_fgamma_baseline(group==2),4.4,0.2,[0,0,0]);
raincloudPlotForPaper(ax,proportion_hfos_baseline(group==1),5,0.2,[1,0,0]);
raincloudPlotForPaper(ax,proportion_hfos_baseline(group==2),5.4,0.2,[0,0,0]);
ylim([0 0.5])
ax.XTick = [0.2,1.2,2.2,3.2,4.2,5.2];
ax.XTickLabel = {'Delta','Theta','Beta','SGamma','FGamma','HFO'};
ylabel('Proportion of power')
xlabel('Frequency band')
title(sprintf('Mean proportion of power per frequency bands - good orgs: %d',onlyGoodOrgs))
% TODO - saveas



cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,kurtosis_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,kurtosis_baseline(group==2),1,0.5,[0,0,0]);

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,skewness_baseline(group==1),0,0.5,[1,0,0]);
raincloudPlotForPaper(ax,skewness_baseline(group==2),1,0.5,[0,0,0]);

%keyboard;

%% Overall changes in power after KCl
%%% overt changes in power/std in 15 minute windows after KCl

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(all_std_15min(group==1,1)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(all_std_15min(group==2,1)),1,0.25,[0,0,0]);

raincloudPlotForPaper(ax,log2(mean(all_std_15min(group==1,2:8),2)),0.5,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(all_std_15min(group==2,2:8),2)),1.5,0.25,[0,0,0]);

ylabel('Log2 Mean StDev 15 minute chunks')
ax.XTick = [0,0.5,1,1.5];
ax.XTickLabel = {'Pre-KCl','Post-KCl','Pre-KCl','Post-KCl'};
xlabel('Group')
title(sprintf('Mean StDev of 15 min chunks pre/post-KCl - goodOrgs=%d',onlyGoodOrgs))
ylim([-8 3])
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3510 0.0565 0.2417 0.7769])
saveas(gcf,sprintf('%s\\MeanRawStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,logscale.png',savepath,onlyGoodOrgs));

ttest2(all_std_15min(group==1,1),mean(all_std_15min(group==1,2:8),2));

%%% Per differentiation


cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(all_std_15min(p1_inds,1)),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(all_std_15min(p1_inds,2:8),2)),0.2,0.1,[1,0,0]);

raincloudPlotForPaper(ax,log2(all_std_15min(p2_inds,1)),0.4,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(all_std_15min(p2_inds,2:8),2)),0.6,0.1,[1,0,0]);

raincloudPlotForPaper(ax,log2(all_std_15min(p3_inds,1)),0.8,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(all_std_15min(p3_inds,2:8),2)),1.0,0.1,[1,0,0]);

raincloudPlotForPaper(ax,log2(all_std_15min(c5_inds,1)),1.2,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean(all_std_15min(c5_inds,2:8),2)),1.4,0.1,[0,0,0]);

raincloudPlotForPaper(ax,log2(all_std_15min(x1_inds,1)),1.6,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean(all_std_15min(x1_inds,2:8),2)),1.8,0.1,[0,0,0]);

ylim([-8 2])
ylabel('Mean Raw StDev 15 minute chunks')
ax.XTick = [0.1,0.5,0.9,1.3,1.7];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
xlabel('Group')
title(sprintf('Mean StDev of 15 min chunks pre/post-KCl - goodOrgs=%d',onlyGoodOrgs))
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259])
saveas(gcf,sprintf('%s\\MeanRawStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,logscale,per-line.png',savepath,onlyGoodOrgs));

%%% Bespoke - only P1 and X1

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(all_std_15min(p1_inds,1)),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(all_std_15min(p1_inds,2:8),2)),0.2,0.1,[1,0,0]);

raincloudPlotForPaper(ax,log2(all_std_15min(x1_inds,1)),0.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean(all_std_15min(x1_inds,2:8),2)),0.6,0.1,[0,0,0]);

ylim([-8 2])
ylabel('Mean Raw StDev 15 minute chunks')
ax.XTick = [0.1,0.5];
ax.XTickLabel = {'P1','X1'};
xlabel('Group')
[~,tmp_p] = ttest2(mean(all_std_15min(p1_inds,2:8),2),mean(all_std_15min(x1_inds,2:8),2)); 
title({sprintf('Log2 Mean StDev of 15 min chunks pre/post-KCl - goodOrgs=%d',onlyGoodOrgs),sprintf('p = %1.2f',tmp_p)})
clear tmp_p
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3870,0.2407,0.1813,0.6259])
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MeanRawStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,logscale,per-line.png',savepath,onlyGoodOrgs));



%% Normalised std now
norm_std_15mins = all_std_15min ./ repmat(all_std_15min(:,1),1,size(all_std_15min,2));
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,norm_std_15mins(group==1,1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,norm_std_15mins(group==2,1),0.5,0.25,[0,0,0]);

raincloudPlotForPaper(ax,mean(norm_std_15mins(group==1,2:8),2),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean(norm_std_15mins(group==2,2:8),2),1.5,0.25,[0,0,0]);

ylim([0 5])
ylabel('Mean StDev 15 minute chunks')
ax.XTick = [0,0.5,1,1.5];
ax.XTickLabel = {'Pre-KCl','Pre-KCl','Post-KCl','Post-KCl'};
xlabel('Group')
title(sprintf('Mean Norm StDev of 15 min chunks pre/post-KCl - goodOrgs=%d',onlyGoodOrgs))
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3510 0.0565 0.2417 0.7769])

saveas(gcf,sprintf('%s\\MeanNormStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,linearscale,ylim[0,5].png',savepath,onlyGoodOrgs));
ylim([0 14])
saveas(gcf,sprintf('%s\\MeanNormStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,linearscale,ylim[0,14].png',savepath,onlyGoodOrgs));


cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(norm_std_15mins(group==1,1)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(norm_std_15mins(group==2,1)),0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(group==1,2:8),2)),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(group==2,2:8),2)),1.5,0.25,[0,0,0]);
ylim([-1 4])
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3302 0.2815 0.2917 0.4102])
xlabel('Group - Pre/Post KCl')
ylabel('Log2 Mean Normalised StDev')
ax.XTickLabel = {'Pre','Post','Pre','Post'};
title('Mean StDev from 15 min windows before and after KCl')
saveas(gcf,sprintf('%s\\MeanNormStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,logscale.png',savepath,onlyGoodOrgs));

%{
%%% Median changes in std in 15 sec windows
cla
ax = gca;
tmp = median(normstds_median(1:60,group==1),1);
raincloudPlotForPaper(ax,tmp(:),1.5,0.25,[1,0,0]);
tmp = median(normstds_median(61:end,group==1),1);
raincloudPlotForPaper(ax,tmp(:),2,0.25,[1,0,0]);
tmp = median(normstds_median(1:60,group==2),1);
raincloudPlotForPaper(ax,tmp(:),2.5,0.25,[0,0,0]);
tmp = median(normstds_median(61:end,group==2),1);
raincloudPlotForPaper(ax,tmp(:),3,0.25,[0,0,0]);
clear tmp
ax.XTick = [1.5,2,2.5,3];
ax.XTickLabel = {'Pre','Post','Pre','Post'}
ylim([0 7])
ylabel('Median normalised StDev')
xlabel('Group - Pre/Post KCl')
title('Median NormStDev from 15 sec windows before and after KCl')
saveas(gcf,sprintf('%s\\MeanNormStDev(15secwindows)Pre-Post-KCl,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));
%}

%%% TODO - median fold-change in the PSD!

%%% CRISPR ones

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,norm_std_15mins(group==1,1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,norm_std_15mins(group==2&crispr_inds,1),0.5,0.25,[0,0,0]);

raincloudPlotForPaper(ax,mean(norm_std_15mins(group==1,2:8),2),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean(norm_std_15mins(group==2&crispr_inds,2:8),2),1.5,0.25,[0,0,0]);

ylim([0 5])
ylabel('Mean StDev 15 minute chunks')
ax.XTick = [0,0.5,1,1.5];
ax.XTickLabel = {'Pre-KCl','Pre-KCl','Post-KCl','Post-KCl'};
xlabel('Group')
title(sprintf('Mean StDev of 15 min chunks pre/post-KCl - goodOrgs=%d',onlyGoodOrgs))
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3510 0.0565 0.2417 0.7769])
saveas(gcf,sprintf('%s\\CRISPR-MeanNormStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,linearscale,ylim[0,5].png',savepath,onlyGoodOrgs));

%%% Per differentiation

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(norm_std_15mins(p1_inds,1)),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(p1_inds,2:8),2)),0.2,0.1,[1,0,0]);

raincloudPlotForPaper(ax,log2(norm_std_15mins(p2_inds,1)),0.4,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(p2_inds,2:8),2)),0.6,0.1,[1,0,0]);

raincloudPlotForPaper(ax,log2(norm_std_15mins(p3_inds,1)),0.8,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(p3_inds,2:8),2)),1.0,0.1,[1,0,0]);

raincloudPlotForPaper(ax,log2(norm_std_15mins(c5_inds,1)),1.2,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(c5_inds,2:8),2)),1.4,0.1,[0,0,0]);

raincloudPlotForPaper(ax,log2(norm_std_15mins(x1_inds,1)),1.6,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(x1_inds,2:8),2)),1.8,0.1,[0,0,0]);

ylim([-1 5])
ylabel('Log2 Mean StDev 15 minute chunks')
ax.XTick = [0.1,0.5,0.9,1.3,1.7];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
xlabel('Group')
title(sprintf('Mean StDev of 15 min chunks pre/post-KCl - goodOrgs=%d',onlyGoodOrgs))
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259]);%[0.3490 0.1935 0.3229 0.6259])
saveas(gcf,sprintf('%s\\MeanNormStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,logscale,per-line,ylim[0,5].png',savepath,onlyGoodOrgs));

%%% Bespoke - P1 and X1 only


cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(norm_std_15mins(p1_inds,1)),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(p1_inds,2:8),2)),0.2,0.1,[1,0,0]);

raincloudPlotForPaper(ax,log2(norm_std_15mins(x1_inds,1)),0.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean(norm_std_15mins(x1_inds,2:8),2)),0.6,0.1,[0,0,0]);

ylim([-1 5])
ylabel('Log2 Mean StDev 15 minute chunks')
ax.XTick = [0.1,0.5];
ax.XTickLabel = {'P1','X1'};
xlabel('Group')
[~,tmp_p] = ttest2(mean(norm_std_15mins(p1_inds,2:8),2),mean(norm_std_15mins(x1_inds,2:8),2)); 
title({sprintf('Log2 Mean StDev of 15 min chunks pre/post-KCl - goodOrgs=%d',onlyGoodOrgs),sprintf('p = %1.2f',tmp_p)})
clear tmp_p
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3870,0.2407,0.1813,0.6259])
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MeanNormStDev(15minwindows)Pre-Post-KCl,goodOrgs-%d,logscale,per-line,ylim[0,5].png',savepath,onlyGoodOrgs));


%% Norm std change in epoch of interest

stdofinterest = normstds_all(61:120,:);
cla
ax = gca;
hold on;

raincloudPlotForPaper(ax,nanmean(stdofinterest(:,group==1),1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,nanmean(stdofinterest(:,group==2),1),1.5,0.25,[0,0,0]);
ylim([0 8])
xlabel('Group')
ylabel('Median Normalised StDev')
ax.XTick = [1 1.5];
ax.XTickLabel = {'Dravet','Control'};
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3380 0.1241 0.2017 0.6333]);
title('Mean NormStd, 15-30min, 15sec win, mean-def')
saveas(gcf,sprintf('%s\\MeanNormStDev(15secwindows;15-30mins),goodOrgs-%d,MeanNormStd.png',savepath,onlyGoodOrgs));

%%% Bespoke - P1 and X1 only
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,nanmean(stdofinterest(:,p1_inds),1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,nanmean(stdofinterest(:,x1_inds),1),1.5,0.25,[0,0,0]);
ylim([0 3])
xlabel('Group')
ylabel('Median Normalised StDev')
ax.XTick = [1 1.5];
ax.XTickLabel = {'P1','X1'};
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3380 0.1241 0.2017 0.6333]);
[~,tmp_p] = ttest2(nanmean(stdofinterest(:,p1_inds),1),nanmean(stdofinterest(:,x1_inds),1)); 
title({'Mean NormStd, 15-30min, 15sec win, mean-def', sprintf('p = %1.2f',tmp_p)})
clear tmp_p
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MeanNormStDev(15secwindows;15-30mins),goodOrgs-%d,MeanNormStd.png',savepath,onlyGoodOrgs));

%%% MEDIAN

stdofinterest = normstds_median(61:121,:);
cla
ax = gca;
hold on;

raincloudPlotForPaper(ax,nanmean(stdofinterest(:,group==1),1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,nanmean(stdofinterest(:,group==2),1),1.5,0.25,[0,0,0]);
ylim([0 8])
xlabel('Group')
ylabel('Median Normalised StDev')
ax.XTick = [1 1.5];
ax.XTickLabel = {'Dravet','Control'};
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3380 0.1241 0.2017 0.6333]);
title('Mean NormStd, 15-30min, 15sec win, median-def')
saveas(gcf,sprintf('%s\\MeanNormStDev(15secwindows;15-30mins),goodOrgs-%d,MedianNormStd.png',savepath,onlyGoodOrgs));

%%% Bespoke - P1 and X1 only
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,nanmean(stdofinterest(:,p1_inds),1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,nanmean(stdofinterest(:,x1_inds),1),1.5,0.25,[0,0,0]);
ylim([0 8])
xlabel('Group')
ylabel('Median Normalised StDev')
ax.XTick = [1 1.5];
ax.XTickLabel = {'P1','X1'};
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3380 0.1241 0.2017 0.6333]);
[~,tmp_p] = ttest2(nanmean(stdofinterest(:,p1_inds),1),nanmean(stdofinterest(:,x1_inds),1)); 
title({'Mean NormStd, 15-30min, 15sec win, mean-def', sprintf('p = %1.2f',tmp_p)})
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MeanNormStDev(15secwindows;15-30mins),goodOrgs-%d,MedianNormStd.png',savepath,onlyGoodOrgs));



%%% Per differentiation


%% PSD Fold change

cla;
ax = gca;
hold on
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(group==1)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(group==2)),0.5,0.25,[0,0,0]);
ylabel('Log2 Fold-change PSD')
ax.XTick = [0, 0.5];
ax.XTickLabel = {'Dravet','Control'};
title('Log2 Fold-change PSD Post/Pre KCl, whole recording, ~[1-50]Hz')
ylim([-1 10])
saveas(gcf,sprintf('%s\\MedianFoldChangePSD[1-50Hz],WholeRecording,goodOrgs-%d,logscale.png',savepath,onlyGoodOrgs));
ranksum(psd_foldchange_lowfreq(group==1),psd_foldchange_lowfreq(group==2))

cla;
ax = gca;
hold on;
ylabel('Log2 Fold-change PSD')
ax.XTick = [0, 0.5];
ax.XTickLabel = {'Dravet','Control'};
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(group==1)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(group==2)),0.5,0.25,[0,0,0]);
title('Log2 Fold-change PSD Post/Pre KCl, 15 mins post-KCl, ~[1-50]Hz')
ylim([-2 8])
saveas(gcf,sprintf('%s\\MedianFoldChangePSD[1-50Hz],15MinsPostKCl,goodOrgs-%d,logscale.pdf',savepath,onlyGoodOrgs));
ranksum(psd_foldchange_15mins(group==1),psd_foldchange_15mins(group==2))

cla;
ax = gca;
hold on;
ylabel('Log2 Fold-change PSD')
ax.XTick = [0, 0.5];
ax.XTickLabel = {'Dravet','Control'};
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(group==1)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(group==2&crispr_inds)),0.5,0.25,[0,0,0]);
title('Log2 Fold-change PSD Post/Pre KCl, 15 mins post-KCl, ~[1-50]Hz')
ylim([-2 8])
saveas(gcf,sprintf('%s\\CRISPR-MedianFoldChangePSD[1-50Hz],15MinsPostKCl,goodOrgs-%d,logscale.png',savepath,onlyGoodOrgs));

%%% Broken down by patient line

cla;
ax = gca;
hold on;
ylabel('Log2 Fold-change PSD')
ax.XTick = [0, 0.5, 1, 1.5, 2];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(p1_inds)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(p2_inds)),0.5,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(p3_inds)),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(c5_inds)),1.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(x1_inds)),2,0.25,[0,0,0]);
title('Log2 Fold-change PSD Post/Pre KCl, whole recording, ~[1-50]Hz')
ylim([-1 10])
saveas(gcf,sprintf('%s\\MedianFoldChangePSD[1-50Hz],WholeRecording,goodOrgs-%d,per-line,logscale.pdf',savepath,onlyGoodOrgs));

cla;
ax = gca;
hold on;
ylabel('Log2 Fold-change PSD')
ax.XTick = [0, 0.5, 1, 1.5, 2];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(p1_inds)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(p2_inds)),0.5,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(p3_inds)),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(c5_inds)),1.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(x1_inds)),2,0.25,[0,0,0]);
title('Log2 Fold-change PSD Post/Pre KCl, 15 mins post-KCl, ~[1-50]Hz')
ylim([-2 8])
saveas(gcf,sprintf('%s\\MedianFoldChangePSD[1-50Hz],15MinsPostKCl,goodOrgs-%d,per-line,logscale.png',savepath,onlyGoodOrgs));


%%% Bespoke - P1 and X1 only
cla;
ax = gca;
hold on;
ylabel('Log2 Fold-change PSD')
ax.XTick = [0, 0.5];
ax.XTickLabel = {'P1','X1'};
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(p1_inds)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_lowfreq(x1_inds)),0.5,0.25,[0,0,0]);
title('Log2 Fold-change PSD Post/Pre KCl, whole recording, ~[1-50]Hz')
ylim([-1 10])
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MedianFoldChangePSD[1-50Hz],WholeRecording,goodOrgs-%d,per-line,logscale.png',savepath,onlyGoodOrgs));

cla;
ax = gca;
hold on;
ylabel('Log2 Fold-change PSD')
ax.XTick = [0, 0.5];
ax.XTickLabel = {'P1','X1'};
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(p1_inds)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(psd_foldchange_15mins(x1_inds)),0.5,0.25,[0,0,0]);
title('Log2 Fold-change PSD Post/Pre KCl, 15 mins post-KCl, ~[1-50]Hz')
ylim([-2 8])
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MedianFoldChangePSD[1-50Hz],15MinsPostKCl,goodOrgs-%d,per-line,logscale.png',savepath,onlyGoodOrgs));


%%% Show all PSDs

figure;
plot(log(psd_all(group==1,:))','Color',[1,0,0,0.75])
hold on
plot(log(psd_all(group==2,:))','Color',[0,0,0,0.75])

plot(mean(log(psd_all(group==1,:)),1),'Color','r','LineWidth',3)
plot(mean(log(psd_all(group==2,:)),1),'Color','k','LineWidth',3)

%% Max Stdev and time to Max Stdev

[maxstd,maxstdind] = max(normstds_all(61:end,:),[],1);

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(nanmean(maxstd(:,group==1),1)),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(nanmean(maxstd(:,group==2),1)),1.5,0.25,[0,0,0]);
ylim([0 8])
xlabel('Group')
ylabel('Log2 Max Mean Normalised StDev')
ax.XTick = [1 1.5];
ax.XTickLabel = {'Dravet','Control'};
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3380 0.1241 0.2017 0.6333]);
title('Max NormStd, 15sec win, mean-def')
saveas(gcf,sprintf('%s\\MaxNormStd(15secwindow),goodOrgs-%d.pdf',savepath,onlyGoodOrgs));

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,nanmean(maxstdind(:,group==1),1)*(15/60)+15,1,0.25,[1,0,0]); % NB - +15 due to KCl time
raincloudPlotForPaper(ax,nanmean(maxstdind(:,group==2),1)*(15/60)+15,1.5,0.25,[0,0,0]);
ylim([0 240])
xlabel('Group')
ylabel('Time of Max Mean Normalised StDev /min')
ax.XTick = [1 1.5];
ax.XTickLabel = {'Dravet','Control'};
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3380 0.1241 0.2017 0.6333]);
[~,tmp_p] = ttest2(nanmean(maxstdind(:,group==1),1)*(15/60)+15,nanmean(maxstdind(:,group==2),1)*(15/60)+15); 
title('Max NormStd Time, 15sec win, mean-def')
saveas(gcf,sprintf('%s\\MaxNormStdTime(15secwindow),goodOrgs-%d.pdf',savepath,onlyGoodOrgs));

%%% bespoke - P1 and X1 only

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(nanmean(maxstd(:,p1_inds),1)),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(nanmean(maxstd(:,x1_inds),1)),1.5,0.25,[0,0,0]);
ylim([0 10])
xlabel('Group')
ylabel('Log2 Max Mean Normalised StDev')
ax.XTick = [1 1.5];
ax.XTickLabel = {'P1','X1'};
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3380 0.1241 0.2017 0.6333]);
[~,tmp_p] = ttest2(nanmean(maxstd(:,p1_inds),1),nanmean(maxstd(:,x1_inds),1)); 
title({'Max NormStd, 15sec win, mean-def',sprintf('p = %1.2f',tmp_p)})
clear tmp_p
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MaxNormStd(15secwindow),goodOrgs-%d.png',savepath,onlyGoodOrgs));

cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,nanmean(maxstdind(:,p1_inds),1)*(15/60)+15,1,0.25,[1,0,0]); % NB - +15 due to KCl time
raincloudPlotForPaper(ax,nanmean(maxstdind(:,x1_inds),1)*(15/60)+15,1.5,0.25,[0,0,0]);
ylim([0 240])
xlabel('Group')
ylabel('Time of Max Mean Normalised StDev /min')
ax.XTick = [1 1.5];
ax.XTickLabel = {'P1','X1'};
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3380 0.1241 0.2017 0.6333]);
[~,tmp_p] = ttest2(nanmean(maxstdind(:,p1_inds),1)*(15/60)+15,nanmean(maxstdind(:,x1_inds),1)*(15/60)+15); 
title({'Max NormStd Time, 15sec win, mean-def',sprintf('p = %1.2f',tmp_p)})
clear tmp_p
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MaxNormStdTime(15secwindow),goodOrgs-%d.png',savepath,onlyGoodOrgs));

%% Time to first std thresh

std_thres = 3;
thres_duration = 2; % number of consecutive windows needing to be above std_thres to count

normstds_thresh = normstds_all(61:end,:)>std_thres;
normstds_thresh = movmean(normstds_thresh,thres_duration,1).*thres_duration; % this should basically count the 1s present in the window of size thres_duration
for org = 1:size(normstds_thresh,2)
    crossing = find(normstds_thresh(:,org)==thres_duration,1,'first');
    if ~isempty(crossing)
        first_crossing(org) = crossing;
    else first_crossing(org) = NaN;
    end
    clear crossing
end
clear org

figure; 
ax = gca;
hold on
raincloudPlotForPaper(ax,first_crossing(group==1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,first_crossing(group==2),1.5,0.25,[0,0,0]);


std_thres = [2:10];
thres_duration = 2; % number of consecutive windows needing to be above std_thres to count
proportions = NaN(length(std_thres),3);

for std = 1:length(std_thres)
    normstds_thresh = normstds_all(61:end,:)>std_thres(std);
    normstds_thresh = movmean(normstds_thresh,thres_duration,1).*thres_duration; % this should basically count the 1s present in the window of size thres_duration

    for org = 1:size(normstds_thresh,2)
        crossing = find(normstds_thresh(:,org)==thres_duration,1,'first');
        if ~isempty(crossing)
            crossing_ind(org) = crossing;
        else
            crossing_ind(org) = NaN;
        end
        clear crossing
    end
    proportions(std,1) = nansum(~isnan(crossing_ind))/length(crossing_ind);
    proportions(std,2) = nansum(~isnan(crossing_ind(group==1)))/length(crossing_ind(group==1));
    proportions(std,3) = nansum(~isnan(crossing_ind(group==2)))/length(crossing_ind(group==2));
    crispr_proportions(std,1) = nansum(~isnan(crossing_ind(group==2&crispr_inds)))/length(crossing_ind(group==2&crispr_inds));
    clear crossing_ind
end
clear org

cla;
ax = gca;
hold on;
plot(std_thres,proportions(:,1),'Color',[0.6,0.6,0.6],'LineWidth',2.5)
plot(std_thres,proportions(:,2),'Color',[1,0,0],'LineWidth',3)
plot(std_thres,proportions(:,3),'Color',[0,0,0],'LineWidth',3)
%plot(std_thres,crispr_proportions(:,1),'Color',[0,0,1],'LineWidth',3)
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3151 0.4296 0.2682 0.3889]);
ylim([0 1])
xlabel('StDev Threshold')
ylabel('Proportion of organoids')
title('Proportion of organoids reaching normstd threshold (mean-defined)')

saveas(gcf,sprintf('%s\\ProportionOrganoids-meanNormStd,goodOrgs-%d.png',savepath,onlyGoodOrgs));


%% Changes in bandpower through time?



%% Changes in bandpower for different normstds (15 second chunks)

%%% First a nice big fancy plot
 
figure;
cla
ax = gca;
hold on;
% Start with baseline quiescence
baseline_activity = normstds_median(1:60,:)' < 1.5 & normstds_median(1:60,:)'>0.5;
for orgo = 1:size(baseline_activity,1)
    delta_baseline(orgo) = nanmean(proportion_delta(orgo,find(baseline_activity(orgo,:))));
end
raincloudPlotForPaper(ax,delta_baseline(group==1),1,0.21,[1,0,0]);
raincloudPlotForPaper(ax,delta_baseline(group==2),1.42,0.21,[0,0,0]);
clear delta_baseline orgo baseline_activity
for std = 2:6
    
    activity_inds = normstds_median' > std;
    for orgo = 1:size(activity_inds,1)
        delta_activity(orgo) = nanmean(proportion_delta(orgo,activity_inds(orgo,:)));
    end
    
    raincloudPlotForPaper(ax,delta_activity(group==1),std,0.21,[1,0,0]);
    raincloudPlotForPaper(ax,delta_activity(group==2),std+0.42,0.21,[0,0,0]);
    clear activity_inds orgo delta_activity
end
ax.XTick = [1.2,2.2,3.2,4.2,5.2,6.2];
ax.XTickLabel = {'Baseline','>2','>3','>4','>5','>6'};
xlabel('Normalised StDev Threshold')
ylabel('Proportional power, delta band')
ylim([0 0.4])
title('Proportion of Delta vs NormStdDev')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1078,0.5000,0.6396,0.3889]);
saveas(gcf,sprintf('%s\\ProportionalPowerPerStd-DeltaBand,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));

%THETA
cla
ax = gca;
hold on;
% Start with baseline quiescence
baseline_activity = normstds_median(1:60,:)' < 1.5 & normstds_median(1:60,:)'>0.5;
for orgo = 1:size(baseline_activity,1)
    delta_baseline(orgo) = nanmean(proportion_theta(orgo,find(baseline_activity(orgo,:))));
end
raincloudPlotForPaper(ax,delta_baseline(group==1),1,0.21,[1,0,0]);
raincloudPlotForPaper(ax,delta_baseline(group==2),1.42,0.21,[0,0,0]);
clear delta_baseline orgo baseline_activity
for std = 2:6
    
    activity_inds = normstds_median' > std;
    for orgo = 1:size(activity_inds,1)
        delta_activity(orgo) = nanmean(proportion_theta(orgo,activity_inds(orgo,:)));
    end
    
    raincloudPlotForPaper(ax,delta_activity(group==1),std,0.21,[1,0,0]);
    raincloudPlotForPaper(ax,delta_activity(group==2),std+0.42,0.21,[0,0,0]);
    clear activity_inds orgo delta_activity
end
ax.XTick = [1.2,2.2,3.2,4.2,5.2,6.2];
ax.XTickLabel = {'Baseline','>2','>3','>4','>5','>6'};
xlabel('Normalised StDev Threshold')
ylabel('Proportional power, theta band')
ylim([0 0.2])
title('Proportion of Theta vs NormStdDev')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1078,0.5000,0.6396,0.3889]);
saveas(gcf,sprintf('%s\\ProportionalPowerPerStd-ThetaBand,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));

%BETA
cla
ax = gca;
hold on;
% Start with baseline quiescence
baseline_activity = normstds_median(1:60,:)' < 1.5 & normstds_median(1:60,:)'>0.5;
for orgo = 1:size(baseline_activity,1)
    delta_baseline(orgo) = nanmean(proportion_beta(orgo,find(baseline_activity(orgo,:))));
end
raincloudPlotForPaper(ax,delta_baseline(group==1),1,0.21,[1,0,0]);
raincloudPlotForPaper(ax,delta_baseline(group==2),1.42,0.21,[0,0,0]);
clear delta_baseline orgo baseline_activity
for std = 2:6
    
    activity_inds = normstds_median' > std;
    for orgo = 1:size(activity_inds,1)
        activity(orgo) = nanmean(proportion_beta(orgo,activity_inds(orgo,:)));
    end
    
    raincloudPlotForPaper(ax,activity(group==1),std,0.21,[1,0,0]);
    raincloudPlotForPaper(ax,activity(group==2),std+0.42,0.21,[0,0,0]);
    clear activity_inds orgo activity
end
ax.XTick = [1.2,2.2,3.2,4.2,5.2,6.2];
ax.XTickLabel = {'Baseline','>2','>3','>4','>5','>6'};
xlabel('Normalised StDev Threshold')
ylabel('Proportional power, beta band')
ylim([0 0.15])
title('Proportion of Beta vs NormStdDev')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1078,0.5000,0.6396,0.3889]);
saveas(gcf,sprintf('%s\\ProportionalPowerPerStd-BetaBand,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));

%SGAMMA
cla
ax = gca;
hold on;
% Start with baseline quiescence
baseline_activity = normstds_median(1:60,:)' < 1.5 & normstds_median(1:60,:)'>0.5;
for orgo = 1:size(baseline_activity,1)
    delta_baseline(orgo) = nanmean(proportion_sgamma(orgo,find(baseline_activity(orgo,:))));
end
raincloudPlotForPaper(ax,delta_baseline(group==1),1,0.21,[1,0,0]);
raincloudPlotForPaper(ax,delta_baseline(group==2),1.42,0.21,[0,0,0]);
clear delta_baseline orgo baseline_activity
for std = 2:6
    
    activity_inds = normstds_median' > std;
    for orgo = 1:size(activity_inds,1)
        activity(orgo) = nanmean(proportion_sgamma(orgo,activity_inds(orgo,:)));
    end
    
    raincloudPlotForPaper(ax,activity(group==1),std,0.21,[1,0,0]);
    raincloudPlotForPaper(ax,activity(group==2),std+0.42,0.21,[0,0,0]);
    clear activity_inds orgo activity
end
ax.XTick = [1.2,2.2,3.2,4.2,5.2,6.2];
ax.XTickLabel = {'Baseline','>2','>3','>4','>5','>6'};
xlabel('Normalised StDev Threshold')
ylabel('Proportional power, sgamma band')
ylim([0.05 0.2])
title('Proportion of Sgamma vs NormStdDev')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1078,0.5000,0.6396,0.3889]);
saveas(gcf,sprintf('%s\\ProportionalPowerPerStd-SlowGammaBand,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));

%FGAMMA
cla
ax = gca;
hold on;
% Start with baseline quiescence
baseline_activity = normstds_median(1:60,:)' < 1.5 & normstds_median(1:60,:)'>0.5;
for orgo = 1:size(baseline_activity,1)
    delta_baseline(orgo) = nanmean(proportion_fgamma(orgo,find(baseline_activity(orgo,:))));
end
raincloudPlotForPaper(ax,delta_baseline(group==1),1,0.21,[1,0,0]);
raincloudPlotForPaper(ax,delta_baseline(group==2),1.42,0.21,[0,0,0]);
clear delta_baseline orgo baseline_activity
for std = 2:6
    
    activity_inds = normstds_median' > std;
    for orgo = 1:size(activity_inds,1)
        activity(orgo) = nanmean(proportion_fgamma(orgo,activity_inds(orgo,:)));
    end
    
    raincloudPlotForPaper(ax,activity(group==1),std,0.21,[1,0,0]);
    raincloudPlotForPaper(ax,activity(group==2),std+0.42,0.21,[0,0,0]);
    clear activity_inds orgo activity
end
ax.XTick = [1.2,2.2,3.2,4.2,5.2,6.2];
ax.XTickLabel = {'Baseline','>2','>3','>4','>5','>6'};
xlabel('Normalised StDev Threshold')
ylabel('Proportional power, fgamma band')
ylim([0.1 0.25])
title('Proportion of Fgamma vs NormStdDev')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1078,0.5000,0.6396,0.3889]);
saveas(gcf,sprintf('%s\\ProportionalPowerPerStd-FastGammaBand,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));

%HFO
cla
ax = gca;
hold on;
% Start with baseline quiescence
baseline_activity = normstds_median(1:60,:)' < 1.5 & normstds_median(1:60,:)'>0.5;
for orgo = 1:size(baseline_activity,1)
    delta_baseline(orgo) = nanmean(proportion_hfos(orgo,find(baseline_activity(orgo,:))));
end
raincloudPlotForPaper(ax,delta_baseline(group==1),1,0.21,[1,0,0]);
raincloudPlotForPaper(ax,delta_baseline(group==2),1.42,0.21,[0,0,0]);
clear delta_baseline orgo baseline_activity
for std = 2:6
    
    activity_inds = normstds_median' > std;
    for orgo = 1:size(activity_inds,1)
        activity(orgo) = nanmean(proportion_hfos(orgo,activity_inds(orgo,:)));
    end
    
    raincloudPlotForPaper(ax,activity(group==1),std,0.21,[1,0,0]);
    raincloudPlotForPaper(ax,activity(group==2),std+0.42,0.21,[0,0,0]);
    clear activity_inds orgo activity
end
ax.XTick = [1.2,2.2,3.2,4.2,5.2,6.2];
ax.XTickLabel = {'Baseline','>2','>3','>4','>5','>6'};
xlabel('Normalised StDev Threshold')
ylabel('Proportional power, HFO band')
ylim([0 0.4])
title('Proportion of HFO vs NormStdDev')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1078,0.5000,0.6396,0.3889]);
saveas(gcf,sprintf('%s\\ProportionalPowerPerStd-HFOBand,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));



%%% Now do them individually...

patent_ea = normstds_median' > 4;
quiescence = normstds_median' < 1.5 & normstds_median'>0.5;
baseline_ea = normstds_median(1:60,:)' < 1.5 & normstds_median(1:60,:)'>0.5;

% Delta first
for orgo = 1:size(patent_ea,1)
    delta_4(orgo) = nanmean(proportion_delta(orgo,patent_ea(orgo,:)));
    delta_qu(orgo) = nanmean(proportion_delta(orgo,quiescence(orgo,61:end)));
    delta_base(orgo) = nanmean(proportion_delta(orgo,find(baseline_ea(orgo,:))));
end

figure;
cla
ax = gca;
hold on;

raincloudPlotForPaper(ax,delta_base(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,delta_base(group==2),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,delta_qu(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,delta_qu(group==2),0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,delta_4(group==1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,delta_4(group==2),1.5,0.25,[0,0,0]);
ax.XTick = [-0.75,0.25,1.25];
ax.XTickLabel = {'Baseline','Post-KCl, Std = 0.5-1.5','Post-KCl, Std > 4'};
ylim([0 0.4])
xlabel('Group and condition')
ylabel('Delta proportional power at different StDevs')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.1943,0.1806,0.3135,0.5741]);
title('Delta proportional power')
saveas(gcf,sprintf('%s\\ProportionalPower-Delta-BaselineVsEA,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));

% Slow gamma
for orgo = 1:size(patent_ea,1)
    sgamma_4(orgo) = nanmean(proportion_sgamma(orgo,patent_ea(orgo,:)));
    sgamma_qu(orgo) = nanmean(proportion_sgamma(orgo,quiescence(orgo,:)));
    sgamma_base(orgo) = nanmean(proportion_sgamma(orgo,find(baseline_ea(orgo,:))));
end

figure;
cla
ax = gca;
hold on;

raincloudPlotForPaper(ax,sgamma_base(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,sgamma_base(group==2),-0.5,0.25,[0,0,0]);

raincloudPlotForPaper(ax,sgamma_qu(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,sgamma_qu(group==2),0.5,0.25,[0,0,0]);

raincloudPlotForPaper(ax,sgamma_4(group==1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,sgamma_4(group==2),1.5,0.25,[0,0,0]);

% Fast first
for orgo = 1:size(patent_ea,1)
    fgamma_4(orgo) = nanmean(proportion_fgamma(orgo,patent_ea(orgo,:)));
    fgamma_qu(orgo) = nanmean(proportion_fgamma(orgo,quiescence(orgo,:)));
    fgamma_base(orgo) = nanmean(proportion_fgamma(orgo,find(baseline_ea(orgo,:))));
end

figure;
cla
ax = gca;
hold on;

raincloudPlotForPaper(ax,fgamma_base(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,fgamma_base(group==2),-0.5,0.25,[0,0,0]);

raincloudPlotForPaper(ax,fgamma_qu(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,fgamma_qu(group==2),0.5,0.25,[0,0,0]);

raincloudPlotForPaper(ax,fgamma_4(group==1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,fgamma_4(group==2),1.5,0.25,[0,0,0]);

% HFO first
for orgo = 1:size(patent_ea,1)
    hfos_4(orgo) = nanmean(proportion_hfos(orgo,patent_ea(orgo,:)));
    hfos_qu(orgo) = nanmean(proportion_hfos(orgo,quiescence(orgo,:)));
    hfos_base(orgo) = nanmean(proportion_hfos(orgo,find(baseline_ea(orgo,:))));
end

figure;
cla
ax = gca;
hold on;

raincloudPlotForPaper(ax,hfos_base(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,hfos_base(group==2),-0.5,0.25,[0,0,0]);

raincloudPlotForPaper(ax,hfos_qu(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,hfos_qu(group==2),0.5,0.25,[0,0,0]);

raincloudPlotForPaper(ax,hfos_4(group==1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,hfos_4(group==2),1.5,0.25,[0,0,0]);

% difference

%% Spectral stuff



patent_ea = normstds_median' > 4;
quiescence = normstds_median' < 1.5 & normstds_median'>0.5;
baseline_ea = normstds_median(1:60,:)' < 1.5 & normstds_median(1:60,:)'>0.5;

%%% KURTOSIS 

mean_kurtosis_prekcl = nanmean(all_kurtosis(:,1:60),2);
mean_kurtosis_postkcl = nanmean(all_kurtosis(:,61:end),2);
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,mean_kurtosis_prekcl(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_prekcl(group==2),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_postkcl(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_postkcl(group==2),0.5,0.25,[0,0,0]);
ylim([1 10])
ax.XTick = [-0.75 0.25];
ax.XTickLabel = {'Pre-KCl','Post-KCl'};
ylabel('Mean spectral kurtosis (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Spectral kurtosis - pre and post-KCl')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1141,0.2046,0.2307,0.5361]);
saveas(gcf,sprintf('%s\\MeanSpectralKurtosis(15secwindows)pre-post-KCl,goodOrgs-%d.png',savepath,onlyGoodOrgs));


%%%% PER LINE

ax = gca;
cla
hold on;
raincloudPlotForPaper(ax,mean_kurtosis_prekcl(p1_inds),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_postkcl(p1_inds),0.2,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_prekcl(p2_inds),0.4,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_postkcl(p2_inds),0.6,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_prekcl(p3_inds),0.8,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_postkcl(p3_inds),1.0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_prekcl(c5_inds),1.2,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_postkcl(c5_inds),1.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_prekcl(x1_inds),1.6,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_kurtosis_postkcl(x1_inds),1.8,0.1,[0,0,0]);
ylim([1 10])
ax.XTick = [0.1,0.5,0.9,1.3,1.7];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
ylabel('Mean spectral kurtosis (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Spectral kurtosis - pre and post-KCl')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259]);%[0.3490 0.1935 0.3229 0.6259])
saveas(gcf,sprintf('%s\\MeanSpectralKurtosis(15secwindows)pre-post-KCl,goodOrgs-%d,per-line,linearscale.png',savepath,onlyGoodOrgs));

%%%% Against STD

for orgo = 1:size(patent_ea,1)
    kurtosis_ea(orgo) = nanmean(all_kurtosis(orgo,patent_ea(orgo,:)));
    kurtosis_qu(orgo) = nanmean(all_kurtosis(orgo,quiescence(orgo,:)));
    kurtosis_base(orgo) = nanmean(all_kurtosis(orgo,find(baseline_ea(orgo,:))));
end
figure;
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,kurtosis_base(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,kurtosis_base(group==2),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,kurtosis_qu(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,kurtosis_qu(group==2),0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,kurtosis_ea(group==1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,kurtosis_ea(group==2),1.5,0.25,[0,0,0]);

%%% SKEWNESS
for orgo = 1:size(patent_ea,1)
    kurtosis_ea(orgo) = nanmean(all_skewness(orgo,patent_ea(orgo,:)));
    kurtosis_qu(orgo) = nanmean(all_skewness(orgo,quiescence(orgo,:)));
    kurtosis_base(orgo) = nanmean(all_skewness(orgo,find(baseline_ea(orgo,:))));
end
figure;
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,kurtosis_base(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,kurtosis_base(group==2),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,kurtosis_qu(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,kurtosis_qu(group==2),0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,kurtosis_ea(group==1),1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,kurtosis_ea(group==2),1.5,0.25,[0,0,0]);

%%% AVERAGE SKEWNESS OVERALL

mean_skewness_prekcl = nanmean(all_skewness(:,1:60),2);
mean_skewness_postkcl = nanmean(all_skewness(:,61:end),2);
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,mean_skewness_prekcl(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_skewness_prekcl(group==2),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,mean_skewness_postkcl(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_skewness_postkcl(group==2),0.5,0.25,[0,0,0]);
ylim([-0.1 2])
ax.XTick = [-0.75 0.25];
ax.XTickLabel = {'Pre-KCl','Post-KCl'};
ylabel('Mean spectral skewness (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Spectral skewness - pre and post-KCl')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1141,0.2046,0.2307,0.5361]);
saveas(gcf,sprintf('%s\\MeanSpectralSkewness(15secwindows)pre-post-KCl,goodOrgs-%d.png',savepath,onlyGoodOrgs));

%%%% PER LINE

ax = gca;
cla
hold on;
raincloudPlotForPaper(ax,mean_skewness_prekcl(p1_inds),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_skewness_postkcl(p1_inds),0.2,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_skewness_prekcl(p2_inds),0.4,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_skewness_postkcl(p2_inds),0.6,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_skewness_prekcl(p3_inds),0.8,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_skewness_postkcl(p3_inds),1.0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_skewness_prekcl(c5_inds),1.2,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_skewness_postkcl(c5_inds),1.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_skewness_prekcl(x1_inds),1.6,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_skewness_postkcl(x1_inds),1.8,0.1,[0,0,0]);
ylim([-0.1 2])
ax.XTick = [0.1,0.5,0.9,1.3,1.7];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
ylabel('Mean spectral skewness (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Spectral skewness - pre and post-KCl')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259]);%[0.3490 0.1935 0.3229 0.6259])
saveas(gcf,sprintf('%s\\MeanSpectralSkewness(15secwindows)pre-post-KCl,goodOrgs-%d,per-line,linearscale.png',savepath,onlyGoodOrgs));

%%% AVERAGE ENERGY OVERALL

mean_energy_prekcl = nanmean(all_energy(:,1:60),2);
mean_energy_postkcl = nanmean(all_energy(:,61:end),2);
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(group==1)),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(group==2)),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(group==1)),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(group==2)),0.5,0.25,[0,0,0]);
ylim([-5 15])
ax.XTick = [-0.75 0.25];
ax.XTickLabel = {'Pre-KCl','Post-KCl'};
ylabel('Log2 Mean signal energy (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Signal energy - pre and post-KCl')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1141,0.2046,0.2307,0.5361]);
saveas(gcf,sprintf('%s\\MeanSpectralEnergy(15secwindows)pre-post-KCl,goodOrgs-%d,logscale.png',savepath,onlyGoodOrgs));

%%% Per line 
ax = gca;
cla
hold on;
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(p1_inds)),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(p1_inds)),0.2,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(p2_inds)),0.4,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(p2_inds)),0.6,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(p3_inds)),0.8,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(p3_inds)),1.0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(c5_inds)),1.2,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(c5_inds)),1.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(x1_inds)),1.6,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(x1_inds)),1.8,0.1,[0,0,0]);
ylim([-5 15])
ylabel('Log2 Mean signal energy (15 second windows) pre and post-KCl');
ax.XTick = [0.1,0.5,0.9,1.3,1.7];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
xlabel('Group/Pre-Post KCl')
title('Signal energy - pre and post-KCl')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259]);%[0.3490 0.1935 0.3229 0.6259])
saveas(gcf,sprintf('%s\\MeanSpectralEnergy(15secwindows)pre-post-KCl,goodOrgs-%d,per-line,logscale.png',savepath,onlyGoodOrgs));

%%% Bespoke - P1 and X1 only
ax = gca;
cla
hold on;
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(p1_inds)),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(p1_inds)),0.2,0.1,[1,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_prekcl(x1_inds)),0.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,log2(mean_energy_postkcl(x1_inds)),0.6,0.1,[0,0,0]);
ylim([-5 15])
ylabel('Log2 Mean signal energy (15 second windows) pre and post-KCl');
ax.XTick = [0.1,0.5];
ax.XTickLabel = {'P1','X1'};
xlabel('Group/Pre-Post KCl')
title('Signal energy - pre and post-KCl')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259]);%[0.3490 0.1935 0.3229 0.6259])
saveas(gcf,sprintf('%s\\Bespoke-P1,X1-Only\\MeanSpectralEnergy(15secwindows)pre-post-KCl,goodOrgs-%d,per-line,logscale.png',savepath,onlyGoodOrgs));


%%% AVERAGE ENTROPY OVERALL

mean_entropy_prekcl = nanmean(all_entropy(:,1:60),2);
mean_entropy_postkcl = nanmean(all_entropy(:,61:end),2);
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,mean_entropy_prekcl(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_entropy_prekcl(group==2),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,mean_entropy_postkcl(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_entropy_postkcl(group==2),0.5,0.25,[0,0,0]);
ylim([0.8 1])
ax.XTick = [-0.75 0.25];
ax.XTickLabel = {'Pre-KCl','Post-KCl'};
ylabel('Mean signal entropy (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Signal entropy - pre and post-KCl')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1141,0.2046,0.2307,0.5361]);
saveas(gcf,sprintf('%s\\MeanSpectralEntropy(15secwindows)pre-post-KCl,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));

%%%% PER LINE

ax = gca;
cla
hold on;
raincloudPlotForPaper(ax,mean_entropy_prekcl(p1_inds),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_entropy_postkcl(p1_inds),0.2,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_entropy_prekcl(p2_inds),0.4,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_entropy_postkcl(p2_inds),0.6,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_entropy_prekcl(p3_inds),0.8,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_entropy_postkcl(p3_inds),1.0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_entropy_prekcl(c5_inds),1.2,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_entropy_postkcl(c5_inds),1.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_entropy_prekcl(x1_inds),1.6,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_entropy_postkcl(x1_inds),1.8,0.1,[0,0,0]);
ylim([0.8 1])
ax.XTick = [0.1,0.5,0.9,1.3,1.7];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
ylabel('Mean signal entropy (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Signal entropy - pre and post-KCl')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259]);%[0.3490 0.1935 0.3229 0.6259])
saveas(gcf,sprintf('%s\\MeanSpectralEntropy(15secwindows)pre-post-KCl,goodOrgs-%d,per-line,linearscale.png',savepath,onlyGoodOrgs));

%%% MEAN SPECTRAL CREST

mean_crest_prekcl = nanmean(all_crest(:,1:60),2);
mean_crest_postkcl = nanmean(all_crest(:,61:end),2);
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,mean_crest_prekcl(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_crest_prekcl(group==2),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,mean_crest_postkcl(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_crest_postkcl(group==2),0.5,0.25,[0,0,0]);
ylim([0 250])
ax.XTick = [-0.75 0.25];
ax.XTickLabel = {'Pre-KCl','Post-KCl'};
ylabel('Mean signal crest (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Signal crest - pre and post-KCl')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1141,0.2046,0.2307,0.5361]);
saveas(gcf,sprintf('%s\\MeanSpectralCrest(15secwindows)pre-post-KCl,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));
%%%% PER LINE

ax = gca;
cla
hold on;
raincloudPlotForPaper(ax,mean_crest_prekcl(p1_inds),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_crest_postkcl(p1_inds),0.2,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_crest_prekcl(p2_inds),0.4,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_crest_postkcl(p2_inds),0.6,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_crest_prekcl(p3_inds),0.8,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_crest_postkcl(p3_inds),1.0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_crest_prekcl(c5_inds),1.2,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_crest_postkcl(c5_inds),1.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_crest_prekcl(x1_inds),1.6,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_crest_postkcl(x1_inds),1.8,0.1,[0,0,0]);
ylim([0 250])
ax.XTick = [0.1,0.5,0.9,1.3,1.7];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
ylabel('Mean signal crest (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Signal crest - pre and post-KCl')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259]);%[0.3490 0.1935 0.3229 0.6259])
saveas(gcf,sprintf('%s\\MeanSpectralCrest(15secwindows)pre-post-KCl,goodOrgs-%d,per-line,linearscale.png',savepath,onlyGoodOrgs));


%%% SPECTRAL FLATNESS

mean_flatness_prekcl = nanmean(all_flatness(:,1:60),2);
mean_flatness_postkcl = nanmean(all_flatness(:,61:end),2);
cla
ax = gca;
hold on;
raincloudPlotForPaper(ax,mean_flatness_prekcl(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_flatness_prekcl(group==2),-0.5,0.25,[0,0,0]);
raincloudPlotForPaper(ax,mean_flatness_postkcl(group==1),0,0.25,[1,0,0]);
raincloudPlotForPaper(ax,mean_flatness_postkcl(group==2),0.5,0.25,[0,0,0]);
ylim([0 0.2])
ax.XTick = [-0.75 0.25];
ax.XTickLabel = {'Pre-KCl','Post-KCl'};
ylabel('Mean spectral flatness (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Spectral flatness - pre and post-KCl')
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1141,0.2046,0.2307,0.5361]);
saveas(gcf,sprintf('%s\\MeanSpectralFlatness(15secwindows)pre-post-KCl,goodOrgs-%d,linearscale.png',savepath,onlyGoodOrgs));

%%%% PER LINE

ax = gca;
cla
hold on;
raincloudPlotForPaper(ax,mean_flatness_prekcl(p1_inds),0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_flatness_postkcl(p1_inds),0.2,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_flatness_prekcl(p2_inds),0.4,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_flatness_postkcl(p2_inds),0.6,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_flatness_prekcl(p3_inds),0.8,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_flatness_postkcl(p3_inds),1.0,0.1,[1,0,0]);
raincloudPlotForPaper(ax,mean_flatness_prekcl(c5_inds),1.2,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_flatness_postkcl(c5_inds),1.4,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_flatness_prekcl(x1_inds),1.6,0.1,[0,0,0]);
raincloudPlotForPaper(ax,mean_flatness_postkcl(x1_inds),1.8,0.1,[0,0,0]);
ylim([0 0.2])
ax.XTick = [0.1,0.5,0.9,1.3,1.7];
ax.XTickLabel = {'P1','P2','P3','C5','X1'};
ylabel('Mean spectral flatness (15 second windows) pre and post-KCl')
xlabel('Group/Pre-Post KCl')
title('Spectral flatness - pre and post-KCl')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2781,0.1917,0.4521,0.6259]);%[0.3490 0.1935 0.3229 0.6259])
saveas(gcf,sprintf('%s\\MeanSpectralFlatness(15secwindows)pre-post-KCl,goodOrgs-%d,per-line,linearscale.png',savepath,onlyGoodOrgs));

%%% DO ENTROPY, SKEWNESS ETC CORRELATE WITH STD?

%%%% averaged over the whole organoid recording
figure;
cla
ax = gca;
hold on;
baseline_inds = (normstds_median(1:60,:)>0.5 & normstds_median(1:60,:)<1.5)';

all_tmp = normstds_median'; % USE THIS AS A SANITY CHECK TEST! 
for org = 1:size(all_energy,1)
    entropy_baseline_15sec_gp1(org) = nanmean(all_kurtosis(org,baseline_inds(org,:)));
    entropy_baseline_15sec_gp2(org) = nanmean(all_kurtosis(org,baseline_inds(org,:)));
end
raincloudPlotForPaper(ax,entropy_baseline_15sec_gp1(group==1),-1,0.25,[1,0,0]);
raincloudPlotForPaper(ax,entropy_baseline_15sec_gp2(group==2),-0.5,0.25,[0,0,0]);

for std = 1:10
    
    tmp_inds = (normstds_median>std-0.5 & normstds_median<std+0.5)';
    for org = 1:size(all_energy,1)
        entropy_std_15sec_gp1(org) = nanmean(all_kurtosis(org,tmp_inds(org,:)));
        entropy_std_15sec_gp2(org) = nanmean(all_kurtosis(org,tmp_inds(org,:)));
    end
    raincloudPlotForPaper(ax,entropy_std_15sec_gp1(group==1),std-1,0.25,[1,0,0]);
    raincloudPlotForPaper(ax,entropy_std_15sec_gp2(group==2),std-0.5,0.25,[0,0,0]);
    
    clear tmp_inds entropy_std_15sec_gp1 entropy_std_15sec_gp2
end

%%%% TODO - for each 15 sec window independently, check it is still
%%%% significant/similar stuff seen

%% Time spent above certain standard deviation

for std = 1:10
    
    tmp_std = normstds_median>std; % CHANGE TO MEDIAN!
    time_above_sec(:,std) = nansum(tmp_std,1)*15;
    
    tmp_std = normstds_median>std-0.5 & normstds_median<std+0.5;
    time_between_sec(:,std) = nansum(tmp_std,1)*15;

    clear tmp_std
end
timegroup = repmat([1:10],34,1);
dravetgroup = repmat(group,1,10);

figure;
cla 
hold on
errorbar([1:10],nanmedian(time_between_sec(group==1,:),1),mad(time_between_sec(group==1,:),[],1),'Color','r','LineWidth',2);
errorbar([1:10],nanmedian(time_between_sec(group==2,:),1),mad(time_between_sec(group==2,:),[],1),'Color','k','LineWidth',2);
title('Median time between StdDevs')
xlabel('Std Dev')
ylabel('Time /sec')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2625 0.2324 0.2917 0.6204]);
saveas(gcf,sprintf('%s\\TimeBetweenStdDevs,goodOrgs-%d,MedianNormStd.png',savepath,onlyGoodOrgs));
anovan(time_between_sec(:),{timegroup(:),dravetgroup(:)},'model','interaction')


figure;
cla 
hold on
errorbar([1:10],nanmedian(time_above_sec(group==1,:),1),mad(time_above_sec(group==1,:),[],1),'Color','r','LineWidth',2);
errorbar([1:10],nanmedian(time_above_sec(group==2,:),1),mad(time_above_sec(group==2,:),[],1),'Color','k','LineWidth',2);
title('Median time above StdDevs')
xlabel('Std Dev')
ylabel('Time /sec')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.2625 0.2324 0.2917 0.6204]);
saveas(gcf,sprintf('%s\\TimeAboveStdDevs,goodOrgs-%d,MedianNormStd.png',savepath,onlyGoodOrgs));
anovan(time_above_sec(:),{timegroup(:),dravetgroup(:)})
