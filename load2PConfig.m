function Config = load2PConfig(DataFiles)


%% Check input arguments, gets UIGet if no datafiles specified
narginchk(0,1);
if ~exist('DataFiles', 'var') || isempty(DataFiles)
    directory = CanalSettings('DataDirectory');
    [DataFiles,p] = uigetfile({'*.sbx;*.tif;*.imgs'}, 'Choose images file(s) to load', directory, 'MultiSelect', 'on');
    if isnumeric(DataFiles)
        Images = []; return
    elseif iscellstr(DataFiles)
        for index = 1:numel(DataFiles)
            DataFiles{index} = fullfile(p,DataFiles{index});
        end
    else
        DataFiles = {fullfile(p,DataFiles)};
    end
elseif ~iscell(DataFiles)
    DataFiles = {DataFiles};
end

%% Load in config information
numFiles = numel(DataFiles);
ext = cell(numFiles, 1);
for index = 1:numFiles
    [~,~,ext{index}] = fileparts(DataFiles{index});
    switch ext{index}
        case '.sbx'
            Config(index) = parseSbxHeader(DataFiles{index});
        case '.tif'
            Config(index) = parseScimHeader(DataFiles{index});  %this is the relevant line for loading our tiff files
            %Scim = scanimage, our acq software.  The header is saved
            %stupidly, so we have to parse the header string to extract
            %experimental info
        case '.imgs'
            Config(index) = parseImgsHeader(DataFiles{index});
        otherwise
            warning('File type %s not recognized...', ext{index});
            Config = [];
    end
end
