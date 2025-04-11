%example organoids

% CTRL
%ORG006 - KCL time is 600
%ORG021

% DRAVET
% ORG015 - P1 and KCl time is 720
%ORG030
%ORG038
    


[num,txt,raw] = xlsread(xlsfile);

%xlsoutputs = 'E:\orgo\OrgoRampingK\OrgoRampingKOutputs.xlsx';

[tmpdir,~,~] = fileparts(xlsfile);
savepath = [tmpdir,'\AnalysisFigures\P1_Only'];

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
    fprintf('\n')
    
    if isnan(raw{org,orgCol})
        fprintf('skipping...\n')
        continue;
    end
    
    if ~raw{org,includeCol} % if not due to be included in the dataset
        fprintf('skipping...\n')
        continue;
    end

    switch raw{org,orgCol}
        case {'ORG010','ORG021'}
            continue;
            colour = [77,171,255]./255;

        case {'ORG006'}
            colour = [77,171,255]./255;
            ylimit = [-1,1];

        case {'ORG015'}
            
            colour = [255,126,121]./255;
            ylimit = [-0.5, 0.5];


        case {'ORG030'}
            
            colour = [255,126,121]./255;
            ylimit = [-0.2, 0.2];

        case {'ORG034'}

            colour = [77,171,255]./255;
            ylimit = [-0.2, 0.2];

        case {'ORG042'}

            colour = [77,171,255]./255;
            ylimit = [-0.4 0.4];
            %? try filtering for high frequencies here?
            keyboard;

        otherwise 
            continue;
            colour = [0,0,0];
    end
    
    group(org) = raw{org,groupCol};
    goodOrg(org) = raw{org,goodOrgCol};
    line{org}    = raw{org,lineCol};
    
    assert(exist(raw{org,saveFileCol},'file')==2)
    load(raw{org,saveFileCol})
    
    
    fs = organoid.Data.Fs;
    kcl_time_sec = raw{org,kclTimeCol};
    
    
    %% TEMPORARY!!!  if .EDR then fix the scaling issue...
    if strcmpi(organoid.Data.meta.fileType,'EDR')
        disp(' Alert: fixing the EDR scaling!')
        rescaling_factor = (organoid.Data.meta.AD./(organoid.Data.meta.ADCMAX*organoid.Data.meta.YCF2));
        fprintf('... Rescaling Y-axis by %3.2f ... ',1/rescaling_factor);
        organoid.Data.data = organoid.Data.data .* rescaling_factor;
    end

    %% Plot the example trace

    xt = [1:length(organoid.Data.data)]./(organoid.Data.Fs*60);
    h = figure('Units','normalized','Position',[0,0.2,0.5,0.40]);

    hold on;
    plot([raw{org,kclTimeCol},raw{org,kclTimeCol}]./60,ylimit,'Color','r','LineStyle',':','LineWidth',2);
    plot(xt,organoid.Data.data,'Color',colour)
    xlim([0 120])
    ylim(ylimit);

  
    title(sprintf('%s: %d - %s - good: %d',raw{org,orgCol},raw{org,kclTimeCol},raw{org,lineCol},raw{org,goodOrgCol}))
   
    saveas(h,sprintf('D:\\Postdoc-IoN\\OrganoidRecordings\\+55mMKCl\\AnalysisFigures\\Examples\\Selected\\%s-%s-%d.pdf',raw{org,orgCol},raw{org,lineCol},raw{org,goodOrgCol}))

    close(h)
    clear end_ind inds
    
    clear fpath fname data organoid normnormcoastline normcoastline normstds 

end

close all % get rid of any unwanted figures that e.g. pentropy may have dragged up


