function [std_thresholds] = findNormStdevCrossings(norm_std,thresholds,consecutive_crossings)


if ~exist('thresholds','var')
    thresholds = [2:8];
end

%% Find first crossing point >x STD

std_thresholds = NaN(length(thresholds),1);

for std_thresh = 1:length(thresholds)
    stds_bool = double(norm_std>thresholds(std_thresh));
    
    if nansum(stds_bool)>0 % threshold crossing found - but was it long enough to count? (i.e. >length of the minimum consecutive train of threshold crossings = > consecutive_crossings)
        
        endinds = strfind(stds_bool',[1 0]);
        
        if isempty(endinds) % In this case, I think the only possibility is that the stds_bool == 1 only at the very end.
            if stds_bool(end)==1
                if any(~stds_bool(end-3:end))
                    std_thresholds(std_thresh) = -1;
                else 
                    std_thresholds(std_thresh) = find(stds_bool,1);
                end
            end
            continue;
        end
        
        stds_boolsum = cumsum(stds_bool);
        endsums = stds_boolsum(endinds);
        stds_bool(endinds+1) = -[endsums(1);diff(endsums)];
        stds_consecutive_thresh = cumsum(stds_bool);

        if ~isempty(find(stds_consecutive_thresh==consecutive_crossings,1))
            std_thresholds(std_thresh) = find(stds_consecutive_thresh==consecutive_crossings,1);
        else 
            std_thresholds(std_thresh) = -1;
        end
        
    else 
        std_thresholds(std_thresh) = -1; % no threshold crossings found
    end 

end