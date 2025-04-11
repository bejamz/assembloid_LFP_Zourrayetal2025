function [dataArray,metadata] = readCFSFile_64b(filepath,params)

% Inputs:
% - params: struct containing metadata about reading in cfs file
%   - .concat = boolean (0/1) of whether to concatenate the individual frames into a single frame of length (fn x 1) where f is number of frames and n is datapoints per frame
%   - .showChannels = boolean of whether to plot channels as it loads them - probably remove this as it only shows it for each frame...?
%   - .manualCheckLoadin = boolean of whether to stop loadin as it progresses for manual check.  Autochanges showChannels etc to 1.
%   - .channelsToLoad = vector of channels of interest to be loaded in, in that order
% Outputs:
% - dataArray = nPoints x nFrames x nCh array of loaded data stream
% - metadata = huge struct with all the metadata for the recording
%   - roughly organised by metadata.frameStuff(i).channelStuff(j)
%
% Grossly edited from JG Colebatch's original read in script to
% - turn into function with inputs
% - remove the manual input steps so can be used in automated readin functions
% - fix bugs in orignal (especially in overwriting of "data" variable when multiple frames or multiple channels read in at once)
%
% Original comments:
%
% cfsRead64c.m 
% an example of reading a cfs file, from 64bit Matlab
% must be able to find CFS64c.dll*
% *equivalent to CED's 64 bit cfs32.dll
% jg Colebatch 
%  11.10.14: made 64 bit compatible version
%

%% Check and sanitise inputs

% Check params is given, if not use defaults
if ~exist('params','var')
    disp('Warning: No params given as input.  Using default values...')
    params = struct;
end
    
% check each field in params in turn
if ~isfield(params,'concat')
    params.concat = 1;
end
if ~isfield(params,'showChannels')
    params.showChannels = 0;
end
if ~isfield(params,'channelsToLoad')
    params.channelsToLoad = []; % all channels
end
if ~isfield(params,'manualCheckLoadin')
    params.manualCheckLoadin = 0;
end

if params.manualCheckLoadin == 1 
    params.showChannels = 1; % this has to be the case if you want to manually check each frame as it is loaded
    % NOT IMPLEMENTED YET!
end


%% Initialise variables
INT1=0;
WRD1=1;
INT2=2;
WRD2=3;
INT4=4;
RL4=5;
RL8=6;
LSTR=7;
EQUALSPACED=0;
MATRIX=1;
SUBSIDIARY=2;
FILEVAR=0;
DSVAR=1;
READ=0;

dataArray = []; % will be size nPoints x nCh x nFrame


%% Check if mex exists (and is compiled)

assert(exist('matcfs64c')==3)

%% Find file of interest / prompt if no file given or no file of interest found
% JGC seemed to care about the matlab release, so let's check we have a version that can handle the incredibly complicated function uigetfile...

if ~exist(filepath,'file') % if the file doesn't exist let's try to load it in manually...

    uiFileWorks = 0;
    s = which('isMATLABReleaseOlderThan'); %introduced in R2020b
    if ~isempty(s)
        verR12 = isMATLABReleaseOlderThan("R12");
        verR2022b = isMATLABReleaseOlderThan("R2022b");
    end
    clear s
    s = which('verLessThan'); % introduced in R2007a
    if ~isempty(s)
        verR12 = verLessThan('matlab','6.0');
        verR2022b = verLessThan('matlab','9.13');
    end
    
    if ~verR12 % strictly speaking impossible as this function was introduced in R2020b... ver R12 was MATLAB 6.0!!!
        caution('...Why are you still running MATLAB <6.0 you degenerate. \n')
    elseif verR12 && ~verR2022b
        caution('...I wrote these functions on MATLAB R2022b, exceptions may occur in older versions.')
        uiFileWorks = 1;
    else 
        uiFileWorks = 1;
    end
    
    clear s verR12 verR2022b

    if uiFileWorks % hurrah your matlab is newer than v6.0...
        [fName, cfsdir]=uigetfile({'*.cfs','CFS files (*.cfs)'}, 'Choose a file to load');
    else 
        %code block basically ripped from original
        % [JGC comment:] note that specifying initial directory only works for Matlab6
        initcfsdir='C:\'; %[JGC comment:] set default
        cfsdir=uigetdir(initcfsdir,'Choose a directory');
        wd=cd;  % save current
        eval(['cd ' '''' cfsdir '''']); % [JGC comment:] otherwise won't work if space in dir name
        [fName, cfsdir]=uigetfile({'*.cfs','CFS files (*.cfs)'}, 'Choose a file to load');
    end

    filepath =  [cfsdir fName];
    clear cfsdir fName
    assert(exist(filepath,'file')) % check we have a valid file now
end

dataTypes=[];
[~,fname,~] = fileparts(filepath);

if length(filepath) > 255 % windows by default has a MAX_PATH of 255 (260 in Windows 10).  This can be changed but I have not tested MEX with longer than normal filepaths
    error('...Filepath length too long! (>255 characters)');
end    

%% Last check of params
% these all must be fields in params for the read-in to work

assert(isfield(params,'showChannels')) 
assert(isfield(params,'manualCheckLoadin'))

%% Open the file!
[fhandle]=matcfs64c('cfsOpenFile',filepath,READ,0); % read only.  Is a mex64 file

if (fhandle < 0)
   error(['...File opening error: ' int2str(fhandle)]); 
end 

% if file opened OK, display basic metadata info
disp(['File information:']); 
byteSize=matcfs64c('cfsFileSize',fhandle);
disp(['File size: ' int2str(byteSize) ' bytes']);
[time,date,comment]=matcfs64c('cfsGetGenInfo',fhandle);   
[channels,fileVars,DSVars,dataSections]=matcfs64c('cfsGetFileInfo',fhandle);
disp(['File ' fname ' created on ' date ' at ' time]);
if  ~(isempty(comment))
    disp(['comment: ' comment]);
end   
disp([int2str(channels) ' channel(s) ']);
disp([int2str(fileVars) ' file variable(s)' ]);
disp([int2str(DSVars) ' data section variable(s)']);
disp([int2str(dataSections) ' data section(s)']);

%disp('paused');
%pause;

%% Read in the file variables in turn and save to metadata

for i=1:fileVars % first let's iterate over the relevant file variables

    [varSize,varType,varUnits,varDesc]=matcfs64c('cfsGetVarDesc',fhandle,i-1,FILEVAR);
    if varType ~= LSTR
        [varValue]=matcfs64c('cfsGetVarVal',fhandle,(i-1),FILEVAR,0,varType);
    else
        [varValue]=matcfs64c('cfsGetVarVal',fhandle,(i-1),FILEVAR,0,varType,varSize); % needed for LSTR
    end 

    disp(' ');
    disp(['FV' int2str(i-1) ':']);
    disp(['Units: ' varUnits]);
    disp(['Description: ' varDesc]);
    switch varType
          case INT1
              dtype='INT1';
          case WRD1
              dtype='WRD1';
          case INT2
              dtype='INT2';
          case WRD2
              dtype='WRD2';
          case INT4
              dtype='INT4';
          case RL4
              dtype='RL4';
          case RL8
              dtype='RL8';
          case LSTR
              dtype='LSTR';
          otherwise
              dtype='unknown';
    end
    disp(['VarType: ' dtype]); 
    if (varType ~=7)
        disp(['Value: ' int2str(varValue)]);
    else
        disp(['Value: ' varValue]); 
    end 

    %% Write file variables into metadata struct
    %%% TO DO

end % for fileVars

%% Read in the data stream itself for all datasections (or a specified subset)
% now for each dataSection or just 1
%if dataSections > 1
%    dsVec=input('frames to read (vector)? ');
%else
%    dsVec=1;
%end
dsVec = params.channelsToLoad(:);
if ~isempty(dsVec) && (length(dsVec)>dataSections || max(dsVec)>dataSections)
    disp('Invalid frames to be read specified.  Using default (all frames)...')
    dsVec = []; % if invalid vector of frames then just load them all in by default
end
if length(dsVec)==0
    dsVec=1:dataSections;
end


for i=1:length(dsVec)

    % Get DS vars from the datafile and put into metadata file
    for j=1:DSVars
        [varSize,varType,varUnits,varDesc]=matcfs64c('cfsGetVarDesc',fhandle,j-1,DSVAR);
        if varType ~= LSTR
            [varValue]=matcfs64c('cfsGetVarVal',fhandle,(j-1),DSVAR,dsVec(i),varType);
        else
            [varValue]=matcfs64c('cfsGetVarVal',fhandle,(j-1),DSVAR,dsVec(i),varType,varSize); % needed for LSTR
        end 

        switch varType
          case INT1
              dtype='INT1';
          case WRD1
              dtype='WRD1';
          case INT2
              dtype='INT2';
          case WRD2
              dtype='WRD2';
          case INT4
              dtype='INT4';
          case RL4
              dtype='RL4';
          case RL8
              dtype='RL8';
          case LSTR
              dtype='LSTR';
          otherwise
              dtype='unknown';
        end
        metadata.cfsMeta.frameMeta(i).DSVars(j).varNum = j-1;
        metadata.cfsMeta.frameMeta(i).DSVars(j).Description = varDesc;
        metadata.cfsMeta.frameMeta(i).DSVars(j).Units = varUnits;
        metadata.cfsMeta.frameMeta(i).DSVars(j).varSize = varSize;
        metadata.cfsMeta.frameMeta(i).DSVars(j).varType = dtype;
        metadata.cfsMeta.frameMeta(i).DSVars(j).Value = varValue;
        clear dtype
    end % for DSVars
    


    [flagSet]=matcfs64c('cfsDSFlags', fhandle, dsVec(i), READ);  % setit = 0 to read
    metadata.cfsMeta.flagSet = flagSet; 
    clear flagSet

    dSbyteSize=matcfs64c('cfsGetDSSize',fhandle,dsVec(i));
    metadata.cfsMeta.datasectionSizeBytes = dSbyteSize;

    % Load in the datastream for each channel at last...
    for j=1:channels
        
        [startOffset,points,yScale,yOffset,xScale,xOffset]=matcfs64c('cfsGetDSChan',fhandle,j-1,dsVec(i));   
        [channelName,yUnits,xUnits,dataType,dataKind,spacing,other]=matcfs64c('cfsGetFileChan',fhandle,j-1);

        if i==1
            dataTypes=[dataTypes dataType];
        end    

      
        % Get metadata for this channel and put into metadata struct
      switch dataType % convert integer representation to string of datatype
          case INT1
              dtype='INT1';
          case WRD1
              dtype='WRD1';
          case INT2
              dtype='INT2';
          case WRD2
              dtype='WRD2';
          case INT4
              dtype='INT4';    
          case RL4
              dtype='RL4';
          case RL8
              dtype='RL8';
          case LSTR
              dtype='LSTR';
          otherwise
              dtype='unknown';
      end
      switch dataKind
          case 0
              dKind='EqSpaced';
          case 1
              dKind='Matrix';
          case 2
              dKind='Subsiduary';
          otherwise
              dKind='unknown';
      end       

      metadata.cfsMeta.frameMeta(i).channel(j).frameNo = i;
      metadata.cfsMeta.frameMeta(i).channel(j).channelName = channelName;
      metadata.cfsMeta.frameMeta(i).channel(j).dataType = dtype;
      metadata.cfsMeta.frameMeta(i).channel(j).xUnits = xUnits;
      metadata.cfsMeta.frameMeta(i).channel(j).yUnits = yUnits;
      metadata.cfsMeta.frameMeta(i).channel(j).dKind = dKind;
      metadata.cfsMeta.frameMeta(i).channel(j).numPoints = points;
      metadata.cfsMeta.frameMeta(i).channel(j).yScale = yScale;
      metadata.cfsMeta.frameMeta(i).channel(j).yOffset = yOffset;
      metadata.cfsMeta.frameMeta(i).channel(j).xScale = xScale;
      metadata.cfsMeta.frameMeta(i).channel(j).xOffset = xOffset;
      metadata.cfsMeta.frameMeta(i).channel(j).samplingFreq = 1./xScale;

      disp(['Channel ' int2str(j-1) ': name is: ' channelName]);
      disp(['yUnits ' yUnits ' xUnits ' xUnits ' datatype ' dtype]);
      disp(['dataKind ' dKind ' spacing ' int2str(spacing) ' other ' int2str(other)]);
      disp(['No of points ' int2str(points)]);
      if xScale ~= 0
          disp(['Sampling frequency: ' num2str(1/xScale)])
      end    

      % Now load the actual points itself
      if points > 0
          metadata.cfsMeta.frameMeta(i).channel(j).dataPresent = 1;
          startPt=0;
          [data]=matcfs64c('cfsGetChanData',fhandle,j-1,dsVec(i),startPt,points,dataType);
          if length(data) ~= points
              disp(['Only ' int2str(length(data)) ' points read!']);
          end    
          data=(data*yScale)+yOffset;
            dataArray(:,j,i) = data; % nPoints x nCh x nFrames

           if params.showChannels=='y'
               xVar=1:length(data);
               xVar=xVar*xScale+xOffset;
               plot(xVar,data);
               xlabel(xUnits);
               title(['Datasection ' int2str(dsVec(i)) ': channel ' int2str(j-1)]);
               disp('paused');
               pause;
           end % if showChan    
            clear data;
      else
          metadata.cfsMeta.frameMeta(i).channel(j).dataPresent = 1;
          %dataArray(:,j,i) = NaN(size(dataArray,1),1); % This only works if the first channel loaded DOES contain points... otherwise we can't increase the size of the dataArray ever
          disp('No points to show');
      end  
     
    end % forloop over channels
end  % forloop over dataSections  

%{
   kspecify=menu('Specific value to read?','yes','no');
   if kspecify==1
    dataS=input('dataSection to read - first is 1? ');  
    chan=input('channel - first = 0? ');
    pointS=input('point to start - first is 0? ');
    pointNum=input('Number to read? ');
    datapt=matcfs64c('cfsGetChanData',fhandle,chan,dataS,pointS,pointNum,dataTypes(chan+1));
    disp(['(Raw) value(s) are: ' num2str(datapt')]);  % not converted with yScale, yOffset
  end    
%}


ret=matcfs64c('cfsCloseFile',fhandle); % close the file

%% Concatenate into single vector if desired

if params.concat

    % The channels should not change name or what they are recording between frames!
    % Therefore we can naively stitch frames together

    concatArray = NaN(size(dataArray,1)*size(dataArray,3),size(dataArray,2)); % nPoints x nCh x nFrames

    for i = 1:size(dataArray,2)
        
        tmp = [];
        for fr = 1:length(dsVec)      
            tmp = cat(1,tmp(:),dataArray(:,i,fr));
        end
        concatArray(:,i) = tmp;
        clear tmp;
        
    end

    dataArray = concatArray;
    clear concatArray

end



end  % end of function 
