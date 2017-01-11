function config = parseScimHeader(TifFile)
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');

% config.version = 1;

%% Check input arguments
narginchk(0,1);
if ~exist('TifFile', 'var') || isempty(TifFile) % Prompt for file selection
    [TifFile, p] = uigetfile({'*.tif'}, 'Select ''tif'' files to load header info from:', directory);
    if isnumeric(TifFile)
        return
    end
    TifFile = fullfile(p, TifFile);
end


%% Set identifying info
config.type = 'scim';
config.FullFilename = TifFile;
[~, config.Filename, ~] = fileparts(TifFile);

%% Identify header information from file

% Load header
metaData = ScanImageTiffReader(TifFile).metadata(); %SI5.2 support

if isempty(metaData)
    
    header = scim_openTif(TifFile);
    config.header = {header};
    
    % Save important information
    if isfield(header,'scanimage') && header.scanimage.VERSION_MAJOR == 5;
        config.Height = header.scanimage.linesPerFrame;
        config.Width = header.scanimage.pixelsPerLine;
        config.Depth = header.scanimage.stackNumSlices;
        config.ZStepSize = header.scanimage.stackZStepSize;
        config.Channels = size(header.scanimage.channelsSave,1);
        config.FrameRate = 1 / header.scanimage.scanFramePeriod;
        config.ZoomFactor = header.scanimage.zoomFactor;
        config.Frames = header.scanimage.acqNumFrames;
    else
        config.Height = header.acq.linesPerFrame;
        config.Width = header.acq.pixelsPerLine;
        config.Depth = header.acq.numberOfZSlices;
        config.ZStepSize = header.acq.zStepSize;
        config.Channels = header.acq.numberOfChannelsSave;
        config.FrameRate = header.acq.frameRate;
        config.ZoomFactor = header.acq.zoomFactor;
        
        % Determine number of frames
        info = imfinfo(TifFile);
        config.Frames = numel(info)/(config.Channels*config.Depth);
    end
    
    
    

    
else %SI5.2
    header.SI = parseSI5Header(metaData);
    config.header = {header};
    
    config.Height = header.SI.hRoiManager.linesPerFrame;
    config.Width = header.SI.hRoiManager.pixelsPerLine;
    config.Depth = numel(header.SI.hStackManager.zs);%numSlices;
    config.ZStepSize = header.SI.hStackManager.stackZStepSize;
    config.Channels = size(header.SI.hChannels.channelSave,1);
    config.FrameRate = 1 / header.SI.hRoiManager.scanFramePeriod;
    config.ZoomFactor = header.SI.hRoiManager.scanZoomFactor;
    
    if config.Depth>1
        try
            config.Frames = header.SI.hFastZ.numVolumes;
            config.framesPerSlice = header.SI.hStackManager.framesPerSlice;
        catch %in case stack from motor
            config.Frames = header.SI.hStackManager.framesPerSlice;
        end
    else
        config.Frames = header.SI.hStackManager.framesPerSlice;
        try %hFastZ may not always be present
            config.Volumes = header.SI.hFastZ.numVolumes;
        end
    end
    
    config.SIversion = header.SI.VERSION_MAJOR; %scanImage Version
    
    
    
end

    %% DEFAULTS
    % config.MotionCorrected = false;
    % config.info = [];
    config.Precision = 'uint16';
    config.DimensionOrder = {'Height','Width','Channels','Depth','Frames'}; % default
    config.Colors = {'green', 'red'};
    config.size = sizeDimensions(config);
