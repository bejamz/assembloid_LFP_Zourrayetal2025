function [] = writeOrganoidStructToExcel(meta_xls,data_xls)

%%% Excel will be structured as follows:
%%% - MetadataOrganoid
%%% - MetadataAnalysis
%%% - SmallBlockAnalysis (Sheet)
%%% - LargeBlockAnalysis (Sheet)

end_time_mins = 120; % as some files are longer than e.g. 2 hours, this specifies explicitly how long we want them to be and will truncate data to fit

[num,txt,raw] = xlsread(meta_xls);

%% Instantiate a new excel if does not exist

if ~exist(data_xls,'file')
    disp('No data XLS file found.  Creating new file...')
end
xlswrite(data_xls,{'ORGANOID ID'},'LargeBlockAnalysis','A1');
xlswrite(data_xls,{'ORGANOID ID'},'SmallBlockAnalysis','A1');

%% Get lengths of segments from the first organoid + write this to MetadataAnalysis
% this should be consistent across all organoids in the list otherwise we can't combine them for stats, obviously
% this will be checked explicitly in the forloop below

%% Get cols from excel file


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
            kclTimeCol = col;
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

%% Write the data itself for each organoid

for org = 3:size(raw,1)
    
    
    if isnan(raw{org,orgCol})
        continue;
    end
    
    if ~raw{org,includeCol} % if not due to be included in the dataset
        continue;
    end
    
    if ~isempty(raw{org,saveFileCol}) && exist(raw{org,saveFileCol},'file')
        fname = raw{org,saveFileCol};
    else
        fpath = raw{org,fpathCol};
        fname_mainfile = raw{org,fnameCol};
        fname = [fpath,'\',raw{org,orgCol},'.mat'];
    end
    
    assert(exist(fname,'file'))
    
    group(org) = raw{org,groupCol};
    goodOrg(org) = raw{org,20};

    load(fname)

    xlswrite(data_xls,raw(org,1),'LargeBlockAnalysis',sprintf('%s1',xlcolumnletter(org-1))); % org-1 because we are starting from position 2 (i.e. B1) and going onwards from there
    
    data_end_index = ceil(end_time_mins/organoid.LargeBlockAnalysis.Meta.LargeBlockLengthMin); % let's get the number of blocks we need here
    tmp_vector = NaN(data_end_index,1);
    end_index = data_end_index;

    row_index_start = 2; % this is the starting point - we begin at row 2 and go downwards
    
    %% LARGE BLOCK - TOTAL POW
    if isfield(organoid.LargeBlockAnalysis,'BandpowerAll')
        if end_index > length(organoid.LargeBlockAnalysis.BandpowerAll) % this checks if the recording is shorter than expected and if so just extracts the bit we have and writes to the tmp_vector, leaving NaNs in the rest
            end_index = length(organoid.LargeBlockAnalysis.BandpowerAll);
        end
        tmp_vector(1:end_index) = organoid.LargeBlockAnalysis.BandpowerAll(1:end_index);
        % assert(length(tmp_vector))==large_block_length);
        xlswrite(data_xls,{'BandpowerAll'},'LargeBlockAnalysis',sprintf('A%d',row_index_start));
        xlswrite(data_xls,tmp_vector(:),'LargeBlockAnalysis',sprintf('%s%d:%s%d',xlcolumnletter(org-1),row_index_start,xlcolumnletter(org-1),row_index_start+length(tmp_vector)-1));
        row_index_start = row_index_start + length(tmp_vector); % the start index for the next data
    end

    %% DELTA 
    if isfield(organoid.LargeBlockAnalysis,'bandpow_delta')
        if end_index > length(organoid.LargeBlockAnalysis.bandpow_delta)
            end_index = length(organoid.LargeBlockAnalysis.bandpow_delta);
        end
        tmp_vector(1:end_index) = organoid.LargeBlockAnalysis.bandpow_delta(1:end_index);
        % assert(length(tmp_vector))==large_block_length);
        xlswrite(data_xls,{'bandpow_delta'},'LargeBlockAnalysis',sprintf('A%d',row_index_start));
        xlswrite(data_xls,tmp_vector(:),'LargeBlockAnalysis',sprintf('%s%d:%s%d',xlcolumnletter(org-1),row_index_start,xlcolumnletter(org-1),row_index_start+length(tmp_vector)-1));
        row_index_start = row_index_start + length(tmp_vector);
    end
    
    %% THETA
    if isfield(organoid.LargeBlockAnalysis,'bandpow_theta')
        if end_index > length(organoid.LargeBlockAnalysis.bandpow_theta)
            end_index = length(organoid.LargeBlockAnalysis.bandpow_theta);
        end
        tmp_vector(1:end_index) = organoid.LargeBlockAnalysis.bandpow_theta(1:end_index);
        % assert(length(tmp_vector))==large_block_length);
        xlswrite(data_xls,{'bandpow_theta'},'LargeBlockAnalysis',sprintf('A%d',row_index_start));
        xlswrite(data_xls,tmp_vector(:),'LargeBlockAnalysis',sprintf('%s%d:%s%d',xlcolumnletter(org-1),row_index_start,xlcolumnletter(org-1),row_index_start+length(tmp_vector)-1));
        row_index_start = row_index_start + length(tmp_vector);
    end

    %% S GAMMA
    if isfield(organoid.LargeBlockAnalysis,'bandpow_sgamma')
        if end_index > length(organoid.LargeBlockAnalysis.bandpow_sgamma)
            end_index = length(organoid.LargeBlockAnalysis.bandpow_sgamma);
        end
        tmp_vector(1:end_index) = organoid.LargeBlockAnalysis.bandpow_sgamma(1:end_index);
        % assert(length(tmp_vector))==large_block_length);
        xlswrite(data_xls,{'bandpow_sgamma'},'LargeBlockAnalysis',sprintf('A%d',row_index_start));
        xlswrite(data_xls,tmp_vector(:),'LargeBlockAnalysis',sprintf('%s%d:%s%d',xlcolumnletter(org-1),row_index_start,xlcolumnletter(org-1),row_index_start+length(tmp_vector)-1));
        row_index_start = row_index_start + length(tmp_vector);
    end

    %% F GAMMA
    if isfield(organoid.LargeBlockAnalysis,'bandpow_fgamma')
        if end_index > length(organoid.LargeBlockAnalysis.bandpow_fgamma)
            end_index = length(organoid.LargeBlockAnalysis.bandpow_fgamma);
        end
        tmp_vector(1:end_index) = organoid.LargeBlockAnalysis.bandpow_fgamma(1:end_index);
        % assert(length(tmp_vector))==large_block_length);
        xlswrite(data_xls,{'bandpow_fgamma'},'LargeBlockAnalysis',sprintf('A%d',row_index_start));
        xlswrite(data_xls,tmp_vector(:),'LargeBlockAnalysis',sprintf('%s%d:%s%d',xlcolumnletter(org-1),row_index_start,xlcolumnletter(org-1),row_index_start+length(tmp_vector)-1));
        row_index_start = row_index_start + length(tmp_vector);
    end

    %% HFO
    if isfield(organoid.LargeBlockAnalysis,'bandpow_hfos')
        if end_index > length(organoid.LargeBlockAnalysis.bandpow_hfos)
            end_index = length(organoid.LargeBlockAnalysis.bandpow_hfos);
        end
        tmp_vector(1:end_index) = organoid.LargeBlockAnalysis.bandpow_hfos(1:end_index);
        % assert(length(tmp_vector))==large_block_length);
        xlswrite(data_xls,{'bandpow_hfos'},'LargeBlockAnalysis',sprintf('A%d',row_index_start));
        xlswrite(data_xls,tmp_vector(:),'LargeBlockAnalysis',sprintf('%s%d:%s%d',xlcolumnletter(org-1),row_index_start,xlcolumnletter(org-1),row_index_start+length(tmp_vector)-1));
        row_index_start = row_index_start + length(tmp_vector);
    end

    %% SMALL BLOCK STD
    row_start_index = 2;
    data_end_index = floor(end_time_mins*60/organoid.SmallBlockAnalysis.Meta.SmallBlockLengthSec); % let's get the number of blocks we need here
    xlswrite(data_xls,raw(org,1),'SmallBlockAnalysis',sprintf('%s1',xlcolumnletter(org-1))); % org-1 because we are starting from position 2 (i.e. B1) and going onwards from there
    if end_index > length(organoid.SmallBlockAnalysis.stdevs_smallblocks)
        end_index = length(organoid.SmallBlockAnalysis.stdevs_smallblocks);
    end
    tmp_vector(1:end_index) = organoid.SmallBlockAnalysis.stdevs_smallblocks(1:end_index);
    % assert(length(tmp_vector))==large_block_length);
    xlswrite(data_xls,{'STDev Small Blocks'},'SmallBlockAnalysis',sprintf('A%d',row_index_start));
    xlswrite(data_xls,tmp_vector(:),'SmallBlockAnalysis',sprintf('%s%d:%s%d',xlcolumnletter(org-1),row_index_start,xlcolumnletter(org-1),row_index_start+length(tmp_vector)-1));
    row_index_start = row_index_start + length(tmp_vector); % the start index for the next data

end


function colLetter = xlcolumnletter(colNumber)

% Excel formats columns using letters.
% This function returns the letter combination that corresponds to a given
% column number.
% Limited to 702 columns
% TAKEN FROM : https://uk.mathworks.com/matlabcentral/answers/54153-dynamic-ranges-using-xlswrite
% NOTE: This requires the function allcomb (see below) - from https://uk.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin

if( colNumber > 26*27 )
    error('XLCOLUMNLETTER: Requested column number is larger than 702. Need to revise method to work with 3 character columns');
else
    % Start with A-Z letters
    atoz        = char(65:90)';
      % Single character columns are first
      singleChar  = cellstr(atoz);
      % Calculate double character columns
      n           = (1:26)';
      indx        = allcomb(n,n);
      doubleChar  = cellstr(atoz(indx));
      % Concatenate
      xlLetters   = [singleChar;doubleChar];
      % Return requested column
      colLetter   = xlLetters{colNumber};
  end

function A = allcomb(varargin)
% FROM https://uk.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin
% ALLCOMB - All combinations
%    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
%    in the arrays A1, A2, ..., and AN. B is P-by-N matrix where P is the product
%    of the number of elements of the N inputs. 
%    This functionality is also known as the Cartesian Product. The
%    arguments can be numerical and/or characters, or they can be cell arrays.
%
%    Examples:
%       allcomb([1 3 5],[-3 8],[0 1]) % numerical input:
%       % -> [ 1  -3   0
%       %      1  -3   1
%       %      1   8   0
%       %        ...
%       %      5  -3   1
%       %      5   8   1 ] ; % a 12-by-3 array
%
%       allcomb('abc','XY') % character arrays
%       % -> [ aX ; aY ; bX ; bY ; cX ; cY] % a 6-by-2 character array
%
%       allcomb('xy',[65 66]) % a combination -> character output
%       % -> ['xA' ; 'xB' ; 'yA' ; 'yB'] % a 4-by-2 character array
%
%       allcomb({'hello','Bye'},{'Joe', 10:12},{99999 []}) % all cell arrays
%       % -> {  'hello'  'Joe'        [99999]
%       %       'hello'  'Joe'             []
%       %       'hello'  [1x3 double] [99999]
%       %       'hello'  [1x3 double]      []
%       %       'Bye'    'Joe'        [99999]
%       %       'Bye'    'Joe'             []
%       %       'Bye'    [1x3 double] [99999]
%       %       'Bye'    [1x3 double]      [] } ; % a 8-by-3 cell array
%
%    ALLCOMB(..., 'matlab') causes the first column to change fastest which
%    is consistent with matlab indexing. Example: 
%      allcomb(1:2,3:4,5:6,'matlab') 
%      % -> [ 1 3 5 ; 1 4 5 ; 1 3 6 ; ... ; 2 4 6 ]
%
%    If one of the N arguments is empty, ALLCOMB returns a 0-by-N empty array.
%    
%    See also NCHOOSEK, PERMS, NDGRID
%         and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)
% Tested in Matlab R2015a and up
% version 4.2 (apr 2018)
% (c) Jos van der Geest
% email: samelinoa@gmail.com
% History
% 1.1 (feb 2006), removed minor bug when entering empty cell arrays;
%     added option to let the first input run fastest (suggestion by JD)
% 1.2 (jan 2010), using ii as an index on the left-hand for the multiple
%     output by NDGRID. Thanks to Jan Simon, for showing this little trick
% 2.0 (dec 2010). Bruno Luong convinced me that an empty input should
% return an empty output.
% 2.1 (feb 2011). A cell as input argument caused the check on the last
%      argument (specifying the order) to crash.
% 2.2 (jan 2012). removed a superfluous line of code (ischar(..))
% 3.0 (may 2012) removed check for doubles so character arrays are accepted
% 4.0 (feb 2014) added support for cell arrays
% 4.1 (feb 2016) fixed error for cell array input with last argument being
%     'matlab'. Thanks to Richard for pointing this out.
% 4.2 (apr 2018) fixed some grammar mistakes in the help and comments


narginchk(1,Inf) ;
NC = nargin ;
% check if we should flip the order
if ischar(varargin{end}) && (strcmpi(varargin{end}, 'matlab') || strcmpi(varargin{end}, 'john'))
    % based on a suggestion by JD on the FEX
    NC = NC-1 ;
    ii = 1:NC ; % now first argument will change fastest
else
    % default: enter arguments backwards, so last one (AN) is changing fastest
    ii = NC:-1:1 ;
end
args = varargin(1:NC) ;
if any(cellfun('isempty', args)) % check for empty inputs
    warning('ALLCOMB:EmptyInput','One of more empty inputs result in an empty output.') ;
    A = zeros(0, NC) ;
elseif NC == 0 % no inputs
    A = zeros(0,0) ; 
elseif NC == 1 % a single input, nothing to combine
    A = args{1}(:) ; 
else
    isCellInput = cellfun(@iscell, args) ;
    if any(isCellInput)
        if ~all(isCellInput)
            error('ALLCOMB:InvalidCellInput', ...
                'For cell input, all arguments should be cell arrays.') ;
        end
        % for cell input, we use to indices to get all combinations
        ix = cellfun(@(c) 1:numel(c), args, 'un', 0) ;
        
        % flip using ii if last column is changing fastest
        [ix{ii}] = ndgrid(ix{ii}) ;
        
        A = cell(numel(ix{1}), NC) ; % pre-allocate the output
        for k = 1:NC
            % combine
            A(:,k) = reshape(args{k}(ix{k}), [], 1) ;
        end
    else
        % non-cell input, assuming all numerical values or strings
        % flip using ii if last column is changing fastest
        [A{ii}] = ndgrid(args{ii}) ;
        % concatenate
        A = reshape(cat(NC+1,A{:}), [], NC) ;
    end
end

