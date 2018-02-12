% The goal of this function is to transform sets of images which depict
% complex, multicellular motility into velocity vector fields which are 
% representative of that motion. To do so, the function loops through 
% several directories to do the following:
%   - Allow the user to define a region of interest and define an axis of 
%   symmetry of a set of multipage tiff stacks
%   - Perform Particle Image Velocimetry (PIV) on the multipage tiff stacks
%   using PIVLab to output a velocity vector field
%   - Transform the velocity vector field into new coordinates defined by
%   the geometry of the region of interest and the axis of symmetry
%   - Analyze this new "primed" velocity vector field for amplitude, 
%   frequency, etc.
% The program assumes *.tif image files.
%
% Inputs:- mainExperimentDirectory (Optional): Directory containing the raw
%            image data. Input 0 if you don't want to use it. Prompts for
%            directory if one isn't given. Push cancel if you don't want to
%            use it.
%        - mainAnalysisDirectory (Optional):   Directory that will/does 
%            contain the analyzed data. Prompts for directory if one isn't 
%            given.
%
% To do: Save the parameters used for analysis in each fish subfolder
%        load manually currated experiment data (e.g. times), plot data
%        Adjust the video/image buttons to autoscale
%        Make variable fields not accept letters
%        Variable lockdown dependent on which analysis is done (e.g.,
%           spatial scale lockdown only after analysis is run)
%        Import button for importing other completed analysis (e.g., PIV
%           data; remember to lockdown template size)

%% Main Function
function analyzeMotility( varargin )

%% Initialize variables

% Initialize directory variables
if( nargin ~= 2 )
    mainExperimentDirectory = uigetdir(pwd, 'Main experiment directory containing image data'); % Main directory containing the subfolders with the images you want to analyze
    mainAnalysisDirectory = uigetdir(pwd, 'Main directory to contain/currently containing analysis'); % Directory to contain/currently containing the analysis from the images (if folder name <analysisFolderName> isn't there, creates it)
else
    mainExperimentDirectory = varargin{1};
    mainAnalysisDirectory = varargin{2};
end
usingExperimentDirectory = ~strcmp(num2str(mainExperimentDirectory),'0') && ~strcmp(mainExperimentDirectory, mainAnalysisDirectory);
if(~usingExperimentDirectory)
    disp('Warning: No experiment directory chosen. Some functions may not be available and may crash the program if used.');
    mainExperimentDirectory = mainAnalysisDirectory;
end
imageFileType = '*.tif'; % Currently redundant (but still necessary), see analysisVariables{1}
currentAnalysesPerformedFileName = 'currentAnalysesPerformed.mat';
rawPIVOutputName = 'rawPIVOutputName'; % WARNING: Don't change this variable name
maskFileOutputName = 'maskVars'; % WARNING: Don't change this variable name
interpolationOutputName = 'processedPIVOutput'; % WARNING: Don't change this variable name
motilityParametersOutputName = 'motilityParameters'; % WARNING: Don't change this variable name
PIVOutputName = 'PIVAnimation.avi'; % WARNING: Don't change this variable name
fourierBoundsDefaults = {'4', '1'};
PIVVideoParams = [0, 1, 100; 0, 1, 100];
nAnalysisCheckboxTypes = 6;

% Initialize GUI variables
startGUI = [1, 1]; % X and Y location of the GUI corner (current units, default is pixels)
screensize = get(groot, 'Screensize'); % Obtain current computer display dimensions
widthGUI = screensize(3) - startGUI(1) - 1; % Define GUI width
heightGUI = screensize(4) - startGUI(2) - 146; % Define GUI height (146 is an empirical number representing my system tray's height)
GUISize = [startGUI(1), startGUI(2), widthGUI, heightGUI]; % Combine parameters for GUI location and dimensions
panelBufferSpacing = 10; % How much spacing is between each panel of logic in the GUI
panelLineWidth = 1;
panelBevelOffset = 2*panelLineWidth + 1;
panelTitleTextSize = 20;
panelTitleHeights = 28;
optionalScrollBarWidth = 100;
experimentVariablesPanelWidthFraction = 0.26; % This variable will multiply 'widthGUI' to determine how wide the variables section of the GUI is; 0.22 seems to be good
analysisPanelPosition = [panelBufferSpacing/widthGUI, panelBufferSpacing/heightGUI, experimentVariablesPanelWidthFraction - panelBufferSpacing/widthGUI, 1/2 - 2*panelBufferSpacing/heightGUI + panelTitleHeights/heightGUI];
analysisTitlePosition = [analysisPanelPosition(1)*widthGUI + panelBevelOffset, analysisPanelPosition(2)*heightGUI + analysisPanelPosition(4)*heightGUI - panelTitleHeights - panelBevelOffset + 2, analysisPanelPosition(3)*widthGUI - 2*panelBevelOffset + 2, panelTitleHeights];
variablesPanelPosition = [analysisPanelPosition(1), analysisPanelPosition(2) + analysisPanelPosition(4) + panelBufferSpacing/heightGUI, analysisPanelPosition(3), 1 - analysisPanelPosition(4) - 3*panelBufferSpacing/heightGUI];
variablesTitlePosition = [variablesPanelPosition(1)*widthGUI + panelBevelOffset, variablesPanelPosition(2)*heightGUI + variablesPanelPosition(4)*heightGUI - panelTitleHeights - panelBevelOffset + 2, variablesPanelPosition(3)*widthGUI - 2*panelBevelOffset + 2, panelTitleHeights];
processingPanelPosition = [analysisPanelPosition(1) + experimentVariablesPanelWidthFraction, analysisPanelPosition(2), (1 - experimentVariablesPanelWidthFraction) - 2*panelBufferSpacing/widthGUI, 1 - 2*panelBufferSpacing/heightGUI];
processingTitlePosition = [processingPanelPosition(1)*widthGUI + panelBevelOffset, processingPanelPosition(2)*heightGUI + processingPanelPosition(4)*heightGUI - panelTitleHeights - panelBevelOffset + 2, processingPanelPosition(3)*widthGUI - 2*panelBevelOffset + 2, panelTitleHeights];
optionalScrollBarPosition = [widthGUI - 2*panelBufferSpacing - optionalScrollBarWidth, processingTitlePosition(2) + 2*panelTitleHeights/3, optionalScrollBarWidth, 1];
widthSubGUI = processingPanelPosition(3)*widthGUI;
heightSubGUI = processingPanelPosition(4)*heightGUI - processingTitlePosition(4);

% Initialize GUI icon size/color variables
textIconHeight = 18;
textBufferSpacing = 4;
inputIconHeight = 5;
answerFieldDropDownWidth = 80;
answerFieldEditWidth = 70;
textSize = 13;
textFGColor = [0, 0, 0];

% Initialize GUI Colors variables
GUIBoxColor = [0.9, 0.925, 0.95];
panelColor = [0.9, 0.9, 0.9];
panelTitleColor = [0.25, 0.25, 0.3];
panelTitleTextColor = [1, 1, 1];

%% Determine experiment directory structures and which analyses, if any, are already done

% Graceful exit if user cancelled one or more directory requests
if(strcmp(num2str(mainAnalysisDirectory),'0'))
    disp('User did not select an analysis directory: Program aborted.');
    return;
end

% Obtain all directory contents, total number of directories
[mainAnalysisDirectoryContents, mainAnalysisSubDirectoryContentsCell, ~] = obtainDirectoryStructure(mainAnalysisDirectory);
if(usingExperimentDirectory)
    [mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainExperimentDirectoryStructuresCorrect] = obtainDirectoryStructure(mainExperimentDirectory);
else
    mainExperimentDirectoryContents = mainAnalysisDirectoryContents;
    mainExperimentSubDirectoryContentsCell = mainAnalysisSubDirectoryContentsCell;
    mainExperimentDirectoryStructuresCorrect = true;
end
nSubDirectories = size(mainExperimentDirectoryContents, 1);

% Verify that image directory structure is correct, abort if not
if(~mainExperimentDirectoryStructuresCorrect)
    disp('Directory structure incorrect: Program aborted. Hint: Main_Directory->Fish_Directory->Vid_Number_Directory->Tiffs');
    return;
end

% Determine if analysis directory already exists, and if not, make it (keep
% in mind that if it exists then it must have the same structure as the 
% image data directories)
determineIfAnalysisFoldersExistsElseCreateThem;

% Remove any subdirectories without tiffs from our list
if(usingExperimentDirectory)
    verifyTiffsInMainDirectoryStructure;
end

% Open or create the main currentAnalysesPerformed.mat file
[currentAnalysisPerformed, analysisVariables] = openOrCreateCurrentAnalysesPerformedFile;

% if(currentAnalysisPerformed.subDirectories == -1) % Should probably fix
%     disp('User did not select an analysis file and didnt want to make one: Program aborted.');
%     return;
% end

analysisToPerform = currentAnalysisPerformed;

%% Create GUI

% Initialize GUI
f = figure('Visible', 'off', 'Position', GUISize, 'Resize', 'off'); % Create figure
set(f, 'name', 'Motility Analysis GUI', 'numbertitle', 'off'); % Rename figure
a = axes; % Define figure axes
set(a, 'Position', [0, 0, 1, 1]); % Stretch the axes over the whole figure
set(a, 'Xlim', [0, widthGUI], 'YLim', [0, heightGUI]); % Switch off autoscaling
set(a, 'XTick', [], 'YTick', []); % Turn off tick marks

% Create background, panels
rectangle('Position', [0, 0, widthGUI, heightGUI], 'Curvature', 0, 'FaceColor', GUIBoxColor, 'Parent', a);
uipanel('Position', variablesPanelPosition, 'Parent', f, 'BackgroundColor', panelColor, 'BorderWidth', panelLineWidth);
uipanel('Position', analysisPanelPosition, 'Parent', f, 'BackgroundColor', panelColor, 'BorderWidth', panelLineWidth);
procPanelParent = uipanel('Position', processingPanelPosition, 'Parent', f, 'BackgroundColor', panelColor, 'BorderWidth', panelLineWidth);

% Create Labels
experimentVariablesTitle  = uicontrol('Parent',f,...
                          'Style','text',...
                          'String','Variables Panel',...
                          'backgroundcolor',panelTitleColor,...
                          'Position',variablesTitlePosition,... % The plus and minus single digits are because of the etched panel
                          'FontName','Gill Sans',...
                          'ForegroundColor',panelTitleTextColor,...
                          'FontSize',panelTitleTextSize); %#ok removes annoying orange warning squiggle under the variable

experimentAnalysisTitle   = uicontrol('Parent',f,...
                          'Style','text',...
                          'String','Analysis Panel',...
                          'backgroundcolor',panelTitleColor,...
                          'Position',analysisTitlePosition,... % The plus and minus single digits are because of the etched panel
                          'FontName','Gill Sans',...
                          'ForegroundColor',panelTitleTextColor,...
                          'FontSize',panelTitleTextSize); %#ok removes annoying orange warning squiggle under the variable

experimentProcessingTitle = uicontrol('Parent',f,...
                          'Style','text',...
                          'String','Image Processing Control Panel',...
                          'backgroundcolor',panelTitleColor,...
                          'Position',processingTitlePosition,... % The plus and minus single digits are because of the etched panel
                          'FontName','Gill Sans',...
                          'ForegroundColor',panelTitleTextColor,...
                          'FontSize',panelTitleTextSize); %#ok removes annoying orange warning squiggle under the variable

% Generate the layout and controls for all of the videos, output the
% boolean array representative of analysis to do/done, data to use, etc.
generateProcessingControlPanelListing;

% Generate the layout for all of experimental variables
generateVariablesPanelListing;

% Generate the layout for the analysis controls
generateAnalysisPanelListing;

% Display GUI
f.Visible = 'on';

%% Auxiliary functions

% obtainDirectoryStructure returns directory structures given a directory
function [mainDirectoryContents, mainSubDirectoryContentsCell, directoryStructuresCorrect] = obtainDirectoryStructure(mainDirectory)

    % Initialize variables
    directoryStructuresCorrect = true;
    
    % Obtain main directory structure
    mainDirectoryContents = dir(mainDirectory); % Obtain all main directory contents
    mainDirectoryContents(~[mainDirectoryContents.isdir]) = []; % Remove non-directories
    mainDirectoryContents(strncmp({mainDirectoryContents.name}, '.', 1)) = []; % Removes . and .. and hidden files
    nSubDirectories = size(mainDirectoryContents, 1);
    mainSubDirectoryContentsCell = cell(1, nSubDirectories);
    
    % Loop through all sub-directory contents to obtain contents
    if(nSubDirectories > 0)
        for i = 1:nSubDirectories
            
            % Obtain main directory structure
            mainSubDirectoryContents = dir(strcat(mainDirectory, filesep, mainDirectoryContents(i).name)); % Obtain all sub-directory contents
            mainSubDirectoryContents(~[mainSubDirectoryContents.isdir]) = []; % Remove non-directories
            mainSubDirectoryContents(strncmp({mainSubDirectoryContents.name}, '.', 1)) = []; % Removes . and .. and hidden files
            mainSubDirectoryContentsCell{i} = mainSubDirectoryContents;
            
            % Check for contents, exit and show error message if empty
            nSubSubDirectories = size(mainSubDirectoryContents, 1);
            if(nSubSubDirectories == 0)
                directoryStructuresCorrect = false;
            end
            
        end
    else
        % No directories, there is obviously a problem
        directoryStructuresCorrect = false;
    end
    
end

% verifyTiffsInMainDirectoryStructure verifies that tiffs are in a 
% directory and removes any directory in which they are not
function verifyTiffsInMainDirectoryStructure
    
    % Make sure directories contain images of the form <imageFileType>, if
    % not, remove the directories from the list
    for i=1:nSubDirectories
        
        nSubSubDirectories = size(mainExperimentSubDirectoryContentsCell{i}, 1);
        directoriesToRemove = [];
        
        % Check for which directories to remove
        for j=1:nSubSubDirectories
            
            subSubDirectoryTiffs = dir(strcat(mainExperimentDirectory,filesep,mainExperimentDirectoryContents(i).name,filesep,mainExperimentSubDirectoryContentsCell{i}(j).name,filesep,imageFileType));
            
            % If directory does not contain tiffs, remove it from the list
            if(isempty(subSubDirectoryTiffs))
                
                % Record directories to remove
                directoriesToRemove = [directoriesToRemove, j]; %#ok since the number of directories should be "small"
                
                % Warn user
                fprintf('Folder "%s->%s" contains no %s''s and will be ignored.\n', mainExperimentDirectoryContents(i).name, mainExperimentSubDirectoryContentsCell{i}(j).name, imageFileType);
                
            end
            
        end
        
        % Actually remove the directories without tiffs
        if(~isempty(directoriesToRemove))
            
            for j=1:size(directoriesToRemove, 2)
            
                mainExperimentSubDirectoryContentsCell{i}(end-j+1) = []; % Removal is done in reverse order as to not mess up numbering of other entries
            
            end
            
        end
        
    end
    
end

% determineIfAnalysisFoldersExistsElseCreateThem name should be obvious
function determineIfAnalysisFoldersExistsElseCreateThem

    % Initialize variables
    directoryAlreadyExists = true;
    
    % Check if the directories are equal
    directoryContentsSameBool = isequal(mainAnalysisDirectoryContents,mainExperimentDirectoryContents);
    
    % Check if the subdirectories are equal
    nMainExperimentSubDirs = size(mainExperimentSubDirectoryContentsCell, 2);
    nMainAnalysisSubDirs = size(mainAnalysisSubDirectoryContentsCell, 2);
    if( (nMainExperimentSubDirs == nMainAnalysisSubDirs)&&(size(mainAnalysisSubDirectoryContentsCell,2)>0) )
        
        for i=1:nMainExperimentSubDirs
            
            if(~isequal({mainExperimentSubDirectoryContentsCell{i}.name},{mainAnalysisSubDirectoryContentsCell{i}.name}))
                
                subDirectoryContentsSameBool = false;
                directoryAlreadyExists = false;
                
            end
            
        end
        
    else
        subDirectoryContentsSameBool = false;
        directoryAlreadyExists = false;
    end
    
    % If analysis directory doesn't perfectly match, make them
    if(~directoryAlreadyExists)
        
        % Match the contents of each directory
        if(~directoryContentsSameBool)
            
            % Check the ith directory
            for i=1:nMainExperimentSubDirs
                
                % Check if the ith directory exists, and if not, make it
                curDirectoryAlreadyExistsInt = exist(strcat(mainAnalysisDirectory,filesep,mainExperimentDirectoryContents(i).name),'dir');
                if(curDirectoryAlreadyExistsInt ~= 7)
                    mkdir(strcat(mainAnalysisDirectory,filesep,mainExperimentDirectoryContents(i).name));
                end
                
            end
            
        end
            
        % Match the contents of each subdirectory
        if( ~subDirectoryContentsSameBool )
            
            % Loop through the directories
            for i=1:nMainExperimentSubDirs
                
                % And check each jth subdirectory
                for j=1:size(mainExperimentSubDirectoryContentsCell{i},1)
                    % Check if the jth subdirectory exists, and if not, make it
                    curSubDirectoryAlreadyExistsInt = exist(strcat(mainAnalysisDirectory,filesep,mainExperimentDirectoryContents(i).name,filesep,mainExperimentSubDirectoryContentsCell{1,i}(j).name),'dir');
                    if(curSubDirectoryAlreadyExistsInt ~= 7)
                        mkdir(strcat(mainAnalysisDirectory,filesep,mainExperimentDirectoryContents(i).name,filesep,mainExperimentSubDirectoryContentsCell{1,i}(j).name));
                    end
                    
                end
                
                
            end
            
        end
        
    end
    
end

% openOrCreateCurrentAnalysesPerformedFile loads the record of which
% analysis has been performed or, if this is done for the first time,
% creates the folder
function [currentAnalysisPerformed, analysisVariables] = openOrCreateCurrentAnalysesPerformedFile
    
    % Obtain file structure
    analysisFile = dir(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName));
    if(isempty(analysisFile))
        createFile = questdlg('Could not locate the analysis file: Create one?');
    end
    
    % If the structure is non-empty, load the file, otherwise make it
    if(~isempty(analysisFile))
        
        currentAnalysisFile = load(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName)); % WARNING: Do not change this variable name without changing the save string below
        currentAnalysisPerformed = currentAnalysisFile.currentAnalysisPerformed;
        analysisFieldNames = fieldnames(currentAnalysisFile);
        if(length(analysisFieldNames) == 1)
            analysisVariables = {'*.tif','32','5','0.1625','1', fourierBoundsDefaults{1}, fourierBoundsDefaults{2}};
        else
            analysisVariables = currentAnalysisFile.analysisVariables;
        end
        if(size(analysisVariables,2) == 5) % Using old variables; update
            analysisVariables{6} = fourierBoundsDefaults{1};
            analysisVariables{7} = fourierBoundsDefaults{2};
        end
        
    elseif(strcmp(createFile,'Yes'))
        
        % Determine how the file should be organized: structure with
        % booleans for each subfolder representing which analysis is done
        nSubDirectories = size(mainExperimentDirectoryContents,1);
        for i=1:nSubDirectories
            currentAnalysisPerformed(i).directory = mainExperimentDirectoryContents(i).name; %#ok since size of structure is small % WARNING: Do not change this variable name without changing the save strings below
            nSubSubDirectories = size(mainExperimentSubDirectoryContentsCell{1, i}, 1);
            currentAnalysisPerformed(i).bools = false(nSubSubDirectories, nAnalysisCheckboxTypes); %#ok since size of structure is small
            currentAnalysisPerformed(i).bools(:, nAnalysisCheckboxTypes) = true; %#ok since size of structure is small
            % Record the order of the bools
            for j=1:nSubSubDirectories
                currentAnalysisPerformedSubDirs(j).subDirectories = mainExperimentSubDirectoryContentsCell{1, i}(j).name; %#ok since size of structure is small
            end
            currentAnalysisPerformed(i).subDirectories = currentAnalysisPerformedSubDirs; %#ok since size of structure is small
        end
        
        % Save this file for future reference, update after any analysis
        analysisVariables = {'*.tif','32','5','0.1625','1', fourierBoundsDefaults{1}, fourierBoundsDefaults{2}};
        save(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName),'currentAnalysisPerformed','analysisVariables'); % WARNING: If currentAnalysisPerformed name is changed, you'll have to manually change this string  IN MANY LOCATIONS!!!
        
    else
        
        currentAnalysisPerformed.subDirectories = -1;
        
    end
    
end

% generateGUIFolderListing builds the display for all the folders and the
% analyses that has been/will be performed. Color indicates if something
% has been done (red no, green yes), checkmark represents what the user
% wants to do.
function generateProcessingControlPanelListing
    
    % Initialize GUI variables
    loopIndex = 1;
    curLinearizedAnalysisIndex = 0;
    textSize = 15;
    textIconHeight = 18;
    textBufferSpacing = 5; % How much spacing is between each line of text
    checkBoxSpacing = 25;
    checkBoxWidth = 20;
    checkboxHeight = 20;
    subFolderWidth = 40;
    subSubFolderWidth = 80 - subFolderWidth - textBufferSpacing;
    nSubSubDirectories = 0;
    subDirTitleColor = [1, 0.6, 0.1];
    subDirTitleTextColor = [1, 1, 1];
    subSubDirTitleColor = [0.5, 0.9, 1];
    subSubDirTitleTextColor = [1, 1, 1];
    notAnalyzedCheckboxColor = [1, 0.8, 0.8];
    analyzedCheckboxColor = [0.8, 1, 0.8];
    playVideoButtonWidth = 55;
    openAnalysisButtonWidth = 80;
    playPIVVideoButtonWidth = 55;
    playPIVSoundButtonWidth = 60;
    loadVarsButtonWidth = 90;
    %g = uipanel('Parent',procPanel);
    
    % Initialize directory structure variables
    for i=1:nSubDirectories
        nSubSubDirectories = nSubSubDirectories + size(mainExperimentSubDirectoryContentsCell{1, i}, 1);
    end
    colWidth = 3*panelBufferSpacing + 4*textBufferSpacing + subFolderWidth + subSubFolderWidth + playVideoButtonWidth + openAnalysisButtonWidth + playPIVVideoButtonWidth + nAnalysisCheckboxTypes*checkBoxSpacing;
    maxRows = floor((heightSubGUI - textBufferSpacing - panelBufferSpacing - checkboxHeight - textIconHeight)/(panelBufferSpacing + textBufferSpacing + 2*textIconHeight));
    maxCols = floor(widthSubGUI/colWidth);
    neededCols = ceil(nSubSubDirectories/maxRows);
    
    if(neededCols > maxCols)
        buttonOverflow = true;
        nRows = maxRows;
        unshownWidth = colWidth*(neededCols - maxCols);
    elseif(nSubSubDirectories/maxRows <= 1)
        buttonOverflow = false;
        nRows = nSubSubDirectories;
        unshownWidth = 0;
    else
        buttonOverflow = false;
        nRows = maxRows;
        unshownWidth = 0;
    end
    
    subProcPanelPosition = [0, 0, widthSubGUI + unshownWidth, processingPanelPosition(4)];
    subProcPanel = uipanel('Position', subProcPanelPosition, 'Parent', procPanelParent, 'BackgroundColor', panelColor, 'BorderWidth', panelLineWidth);
    
%    nCols = ceil(((textIconHeight + textBufferSpacing)*(nSubDirectories + 2*nSubSubDirectories) + textBufferSpacing)/(heightSubGUI - 2*textIconHeight - 3*textBufferSpacing)); % Minus 2 is for the header and the titles
%    nRows = ceil((nSubDirectories + 2*nSubSubDirectories)/nCols);
    checkBoxHandleArray = gobjects(nSubSubDirectories, nAnalysisCheckboxTypes);
    buttonHandleArray = gobjects(nSubSubDirectories, 6);
    
    % Enumerate each folder and which analyses have been performed/will be
    % performed
    for i=1:nSubDirectories
        
        % Determine how many sub folders are in the current directory
        curSubSubSize = size(mainExperimentSubDirectoryContentsCell{1,i},1);
        
        % Determine where to put the current subfolder label
        curRow = mod(loopIndex-1,nRows) + 1;
        curCol = ceil(loopIndex/nRows);
        
        % Create the current subfolder label
        
        curSubFolderX = panelBufferSpacing + (curCol - 1)*(widthSubGUI + unshownWidth)/neededCols;
        curSubFolderY = heightSubGUI - 2*textIconHeight - 2*panelBufferSpacing - textBufferSpacing - checkboxHeight - (curRow - 1)*(panelBufferSpacing + 2*textIconHeight + textBufferSpacing);
        
        uicontrol('Parent',subProcPanel,...
            'Style','text',...
            'String',mainExperimentDirectoryContents(i).name,...
            'backgroundcolor',subDirTitleColor,...
            'Position',[curSubFolderX, curSubFolderY, subFolderWidth, textIconHeight],... % The plus and minus 1 are because of the etched panel
            'FontName','Gill Sans',...
            'ForegroundColor',subDirTitleTextColor,...
            'FontSize',textSize);
        
        % Determine which analyses to perform (inverse of already done)
        analysisToPerform(i).bools(:,1:end -1) = ~analysisToPerform(i).bools(:,1:end -1); % Flip everything except the 'Use' section, 1 used to mean done, so we flip it to 0 so we don't perform analysis on it, except the 'Use' bool
        
        % Loop through the current folders subdirectories
        for j=1:curSubSubSize
            
            curLinearizedAnalysisIndex = curLinearizedAnalysisIndex + 1;
            curRow = mod(loopIndex-1,nRows) + 1;
            curCol = ceil(loopIndex/nRows);
            curDividingPointStart = (curCol - 1)*(widthSubGUI + unshownWidth)/neededCols;
            curDividingPointEnd = curCol*(widthSubGUI + unshownWidth)/neededCols;
            curXStart = panelBufferSpacing + subFolderWidth + curDividingPointStart;
            curXEnd = curDividingPointEnd - panelBufferSpacing;
            curY = heightSubGUI - 2*textIconHeight - 2*panelBufferSpacing - textBufferSpacing - checkboxHeight - (curRow - 1)*(panelBufferSpacing + 2*textIconHeight + textBufferSpacing);
            curPosition = [curXStart, curY, subSubFolderWidth, textIconHeight];
            
            % Create the current subfolder label
            uicontrol('Parent',subProcPanel,...
                'Style','text',...
                'String',mainExperimentSubDirectoryContentsCell{i}(j).name,...
                'backgroundcolor',subSubDirTitleColor,...
                'Position',curPosition,...
                'FontName','Gill Sans',...
                'ForegroundColor',subSubDirTitleTextColor,...
                'FontSize',textSize);
            
            % Create play video button
            buttonHandleArray(loopIndex, 1) = uicontrol('Parent',subProcPanel,...
                'Style','pushbutton',...
                'String','Play Video',...
                'Position',[curPosition(1) + subSubFolderWidth + textBufferSpacing, curPosition(2), playVideoButtonWidth, textIconHeight],...
                'Callback',{@playVideoButton_Callback, i, j});
%             if( currentAnalysisPerformed(i).bools(j, 4) )
%                 set(buttonHandleArray(loopIndex, 1), 'Enable', 'on');
%             else
%                 set(buttonHandleArray(loopIndex, 1), 'Enable', 'off');
%             end
            
            % Create open images button
            buttonHandleArray(loopIndex, 2) = uicontrol('Parent',subProcPanel,...
                'Style','pushbutton',...
                'String','Analysis Images',...
                'Position',[curPosition(1) + subSubFolderWidth + playVideoButtonWidth + 2*textBufferSpacing, curPosition(2), openAnalysisButtonWidth, textIconHeight],...
                'Callback',{@openAnalysisImagesButton_Callback, i, j});
            if( currentAnalysisPerformed(i).bools(j, 4) )
                set(buttonHandleArray(loopIndex, 2), 'Enable', 'on');
            else
                set(buttonHandleArray(loopIndex, 2), 'Enable', 'off');
            end
            
            % Create open piv video button
            buttonHandleArray(loopIndex, 3) = uicontrol('Parent',subProcPanel,...
                'Style','pushbutton',...
                'String','Play PIV Video',...
                'Position',[curPosition(1) + subSubFolderWidth + playVideoButtonWidth + openAnalysisButtonWidth + 3*textBufferSpacing, curPosition(2), playPIVVideoButtonWidth, textIconHeight],...
                'Callback',{@playPIVVideoButton_Callback, i, j});
            if( currentAnalysisPerformed(i).bools(j, 4) )
                set(buttonHandleArray(loopIndex, 3), 'Enable', 'on');
            else
                set(buttonHandleArray(loopIndex, 3), 'Enable', 'off');
            end
            
            % Create play sound button
            buttonHandleArray(loopIndex, 4) = uicontrol('Parent',subProcPanel,...
                'Style','togglebutton',...
                'String','Play Sound',...
                'Position',[curPosition(1) + subSubFolderWidth + textBufferSpacing, curPosition(2) - (textIconHeight + textBufferSpacing), playPIVSoundButtonWidth, textIconHeight],...
                'Callback',{@playPIVSoundButton_Callback, i, j});
            if( currentAnalysisPerformed(i).bools(j, 4) )
                set(buttonHandleArray(loopIndex, 4), 'Enable', 'on');
            else
                set(buttonHandleArray(loopIndex, 4), 'Enable', 'off');
            end
            
            % Create load into workspace button
            buttonHandleArray(loopIndex, 5) = uicontrol('Parent',subProcPanel,...
                'Style','pushbutton',...
                'String','Vars->Workspace',...
                'Position',[curPosition(1) + subSubFolderWidth + 2*textBufferSpacing + playPIVSoundButtonWidth, curPosition(2) - (textIconHeight + textBufferSpacing), loadVarsButtonWidth, textIconHeight],...
                'Callback',{@loadVarsButton_Callback, i, j});
            if( currentAnalysisPerformed(i).bools(j, 4) )
                set(buttonHandleArray(loopIndex, 5), 'Enable', 'on');
            else
                set(buttonHandleArray(loopIndex, 5), 'Enable', 'off');
            end
            
            % Create all of the checkboxes
            for k=1:nAnalysisCheckboxTypes
                
                currentlyDone = ~analysisToPerform(i).bools(j, k);
                if( k ~= nAnalysisCheckboxTypes )
                    currentColor = currentlyDone*analyzedCheckboxColor + ~currentlyDone*notAnalyzedCheckboxColor; % Simple way of changing colors depending on whether or not analysis is done
                else
                    currentColor = ~currentlyDone*analyzedCheckboxColor + currentlyDone*notAnalyzedCheckboxColor; % Remember to flip the color of the 'Use' section
                end
                % Create the current subfolder label
                checkBoxHandleArray(curLinearizedAnalysisIndex, k) = uicontrol('Parent', subProcPanel,...
                    'Style', 'checkbox',...
                    'BackgroundColor', currentColor,...
                    'Position', [curXEnd - (nAnalysisCheckboxTypes - k + 1)*checkBoxSpacing - 4, curY, checkBoxWidth, checkboxHeight],...
                    'Value',~currentlyDone,...
                    'Callback', {@checkBox_Callback, curLinearizedAnalysisIndex, i, j, k});
                
            end
            
            % Since each subsubdir gets two rows, add another to the index
            loopIndex = loopIndex + 1;
            
        end
        
    end
    
    % Build the select all checkboxes
    for i=1:nAnalysisCheckboxTypes
        
        % Determine if it should be initially checked or not
        currentCheckedCheck = zeros(1, nSubDirectories);
        for j=1:nSubDirectories
            currentCheckedCheck(j) = min(analysisToPerform(j).bools(:, i));
        end
        
        % Create label for subfolder
        textSizeModified = -4;
        curLabelPosition = [(widthSubGUI + unshownWidth)/neededCols - panelBufferSpacing - (nAnalysisCheckboxTypes - i + 1)*checkBoxSpacing - 4, heightSubGUI - textIconHeight - panelBufferSpacing, 0, checkboxHeight];
        
        switch i
            
            case 1
                % Create the current subfolder label
                curLabelWidth = [-15, 0, 25, 0];
                uicontrol('Parent',subProcPanel,...
                    'Style','text',...
                    'String','PIV',...
                    'backgroundcolor',panelColor,...
                    'Position',curLabelPosition + curLabelWidth,...
                    'FontName','Gill Sans',...
                    'FontSize',textSize + textSizeModified);
                
            case 2
                
                % Create the current subfolder label
                curLabelWidth = [-22, 0, 45, 0];
                uicontrol('Parent',subProcPanel,...
                    'Style','text',...
                    'String','Outline',...
                    'backgroundcolor',panelColor,...
                    'Position',curLabelPosition + curLabelWidth,...
                    'FontName','Gill Sans',...
                    'FontSize',textSize + textSizeModified);
                
            case 3
                % Create the current subfolder label
                curLabelWidth = [-12, 0, 35, 0];
                uicontrol('Parent',subProcPanel,...
                    'Style','text',...
                    'String','Interp',...
                    'backgroundcolor',panelColor,...
                    'Position',curLabelPosition + curLabelWidth,...
                    'FontName','Gill Sans',...
                    'FontSize',textSize + textSizeModified);
                
            case 4
                % Create the current subfolder label
                curLabelWidth = [-8, 0, 40, 0];
                uicontrol('Parent',subProcPanel,...
                    'Style','text',...
                    'String','Analyze',...
                    'backgroundcolor',panelColor,...
                    'Position',curLabelPosition + curLabelWidth,...
                    'FontName','Gill Sans',...
                    'FontSize',textSize + textSizeModified);
                
            case 5
                % Create the current subfolder label
                curLabelWidth = [0, 0, 32, 0];
                uicontrol('Parent',subProcPanel,...
                    'Style','text',...
                    'String','Video',...
                    'backgroundcolor',panelColor,...
                    'Position',curLabelPosition + curLabelWidth,...
                    'FontName','Gill Sans',...
                    'FontSize',textSize + textSizeModified);
                
            case 6
                % Create the current subfolder label
                curLabelWidth = [0, 0, 25, 0];
                uicontrol('Parent',subProcPanel,...
                    'Style','text',...
                    'String','Use',...
                    'backgroundcolor',panelColor,...
                    'Position',curLabelPosition + curLabelWidth,...
                    'FontName','Gill Sans',...
                    'FontSize',textSize + textSizeModified);
                
            otherwise
                % Create the current subfolder label
                curLabelWidth = [-15, 0, 45, 0];
                uicontrol('Parent',subProcPanel,...
                    'Style','text',...
                    'String','Unknown',...
                    'backgroundcolor',panelColor,...
                    'Position',curLabelPosition + curLabelWidth,...
                    'FontName','Gill Sans',...
                    'FontSize',textSize + textSizeModified);
                
        end
        
        % Create the current subfolder label (only checked if all entries
        % down the rows are checked)
        curCheckboxPosition = [curLabelPosition(1), curLabelPosition(2) - textBufferSpacing - textIconHeight, checkBoxWidth, checkboxHeight];
        selectAllHandle = uicontrol('Parent', subProcPanel,...
            'Style', 'checkbox',...
            'Position', curCheckboxPosition,...
            'backgroundcolor',panelColor,...
            'Value',min(currentCheckedCheck),...
            'Callback', {@selectAllCheckBox_Callback,i});
        
        if( i == nAnalysisCheckboxTypes )
            set(selectAllHandle, 'Enable', 'off');
        end
        
    end
    
    if(buttonOverflow)
        
        panelSlider = uicontrol('Parent', f,...
            'Style', 'slider',...
            'Position', optionalScrollBarPosition,...
            'backgroundcolor',panelColor,...
            'Value',0,...
            'Callback', {@panelSlider_Callback}); %#ok
    end
        
    % Callback functions
    % Checkbox functionality
    function checkBox_Callback(~, ~, curLinArrayIndex, ii, jj, kk)
        
        % Change boolean value
        analysisToPerform(ii).bools(jj, kk) = ~analysisToPerform(ii).bools(jj, kk);
        
        % Update the color/file if we are toggling the 'Use' checkbox
        if(kk == nAnalysisCheckboxTypes)
            
            % Change the color of the checkbox
            curColor = analysisToPerform(ii).bools(jj, kk)*analyzedCheckboxColor + ~analysisToPerform(ii).bools(jj, kk)*notAnalyzedCheckboxColor;
            set(checkBoxHandleArray(curLinArrayIndex, kk), 'BackgroundColor', curColor);
            
            % Update the file
            currentAnalysisPerformed(ii).bools(jj, kk) = analysisToPerform(ii).bools(jj, kk);
            % Save this file for future reference, update after any analysis
            save(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName),'currentAnalysisPerformed','analysisVariables'); % WARNING: If currentAnalysisPerformed name is changed, you'll have to manually change this string IN MANY LOCATIONS!!!
            
        end
        
    end
    
    % Select all checkbox functionality
    function selectAllCheckBox_Callback( hObject, ~, kk)
        
        % Check if boxes in row are mixed on and off
        uniqueEntries = [];
        for ii=1:nSubDirectories
            uniqueEntries = unique([unique(analysisToPerform(ii).bools(:, kk)); uniqueEntries]);
        end
        
        % If some are on and others off, turn all on, else toggle
        if(size(uniqueEntries, 1) > 1)
            
            % Turn analysis to perform on
            for aTPIndex=1:size(analysisToPerform, 2)
                analysisToPerform(aTPIndex).bools(:, kk) = true;
            end
            
            % Set current value
            set(hObject, 'Value', 1);
            
            % Set all checkboxes in column to true
            for ii=1:nSubSubDirectories
                set(checkBoxHandleArray(ii, kk), 'Value', 1)
            end
            
        else
            
            % Get the current value
            curBoolValue = get(hObject, 'Value');
            
            % Set analysisToPerform to opposite its current value
            for aTPIndex=1:size(analysisToPerform, 2)
                analysisToPerform(aTPIndex).bools(:, kk) = ~analysisToPerform(aTPIndex).bools(:, kk);
            end
            
            % Set all checkboxes in column to true
            for ii=1:nSubSubDirectories
                set(checkBoxHandleArray(ii, kk), 'Value', curBoolValue)
            end
            
        end
        
    end

    function playVideoButton_Callback(~, ~, ii, jj)

        % Define the current folder, set variables
        curExpDir = strcat(mainExperimentDirectory, filesep, mainExperimentDirectoryContents(ii).name, filesep, mainExperimentSubDirectoryContentsCell{1, ii}(jj).name);
        
        % Obtain current size of combined movies
        dirContents = dir(strcat(curExpDir,filesep,'*.tif'));
        totalFileSize = sum([dirContents.bytes]);
        timeReducedFileSize = totalFileSize*(PIVVideoParams(1,3) - PIVVideoParams(1,1))/(PIVVideoParams(1,2)*100);
        curFileSize = timeReducedFileSize*(PIVVideoParams(2,3) - PIVVideoParams(2,1))/(PIVVideoParams(2,2)*100);
        curFileSize = curFileSize/(1e9);
        
        % Warn user that file may be too large
        questAns = questdlg(sprintf('This function will load all parts of the movie specified. This may result in loading large volumes of data into RAM and possible crashing. With current options, this will load roughly %.3g GB into RAM. Continue?\n\n Note: Settings can be changed in the lower left panel under various video settings.',curFileSize),'File Size Warning','Yes','No','Yes');
        
        if(strcmp(questAns,'Yes'))
            % Create Video
            createNonPIVMovie(curExpDir, PIVVideoParams);
        end
        
    end

    function openAnalysisImagesButton_Callback(~, ~, ii, jj)
        
        % Define the current folder
        curDir = strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(ii).name, filesep, mainExperimentSubDirectoryContentsCell{1, ii}(jj).name);
        
        % Load the figure
        openfig(strcat(curDir, filesep, 'Figures_Current.fig'));
        
    end

    function playPIVVideoButton_Callback(~, ~, ii, jj)
        
        % Define the current folder, set variables
        curExpDir = strcat(mainExperimentDirectory, filesep, mainExperimentDirectoryContents(ii).name, filesep, mainExperimentSubDirectoryContentsCell{1, ii}(jj).name);
        curAnDir = strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(ii).name, filesep, mainExperimentSubDirectoryContentsCell{1, ii}(jj).name);
        questAns = 'N/A';
        
        % If the file has not been created yet, create it
        if(~currentAnalysisPerformed(ii).bools(jj, 5))
            questAns = questdlg('Video has not been made yet: create now?','No video exists','Yes','No','Yes');
        end
        
        % If the file exists (N/A), load it, otherwise create it or don't
        if(strcmp(questAns,'N/A'))
            
            implay(strcat(curAnDir,filesep,PIVOutputName))
            
        elseif(strcmp(questAns,'Yes'))
            
            % Close GUI and everything else
            close all;
            
            % Create Video
            createPIVMovie(curAnDir, curExpDir, analysisVariables, PIVVideoParams, PIVOutputName, interpolationOutputName);
            
            % Update currentAnalysisPerformed
            currentAnalysisPerformed(ii).bools(jj,5) = true;
            
            % Save currentAnalysisPerformed.mat
            save(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName),'currentAnalysisPerformed','analysisVariables'); % WARNING: If currentAnalysisPerformed name is changed, you'll have to manually change this string IN MANY LOCATIONS!!!
            
            % Reopen this program
            analyzeMotility(mainExperimentDirectory, mainAnalysisDirectory);
            
        end
        
    end
    
    function playPIVSoundButton_Callback(hObject, ~, ii, jj)
        
        % Get button state
        isToggleDown = get(hObject, 'Value');
        set(hObject, 'String', 'Loading...');
        pause(0.1);
        
        % Define the current folder
        curDir = strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(ii).name, filesep, mainExperimentSubDirectoryContentsCell{1, ii}(jj).name);

        if(logical(isToggleDown))
            
            % Load data
            gutMeshVelsPCoordsStruct = load(strcat(curDir, filesep, 'processedPIVOutput_Current.mat'),'gutMeshVelsPCoords');
            gutMeshVelsPCoords = gutMeshVelsPCoordsStruct.gutMeshVelsPCoords;
            
            % Play sound if toggle was pressed down, stop if released
            totalTimeFraction = 1;
            fractionOfTimeStart = size(gutMeshVelsPCoords,4);
            markerNumStart = 1;
            markerNumEnd = size(gutMeshVelsPCoords,2);
            samplingRate = 44100;
            
            % Transform into sound with position mapped onto frequency
            theSound = playMotilityAsSound(gutMeshVelsPCoords, totalTimeFraction, fractionOfTimeStart, markerNumStart, markerNumEnd);
            
            % Set string
            set(hObject, 'String', 'Stop Sound');
            
            % Play the sound
            sound(theSound,samplingRate)
            
        else
            
            % Stop the sound
            clear sound
            
            % Set string
            set(hObject, 'String', 'Play Sound');
            
        end
        
    end
    
    function loadVarsButton_Callback(~, ~, ii, jj)
        
        % Define the current folder
        curDir = strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(ii).name, filesep, mainExperimentSubDirectoryContentsCell{1, ii}(jj).name);
        
        % Load the variables into the workspace
        assignin('base', 'curDir', curDir);
        evalin('base', 'load(strcat(curDir, filesep, ''motilityParameters_Current.mat''))');
        evalin('base', 'load(strcat(curDir, filesep, ''processedPIVOutput_Current.mat''))');
        
    end
    
    function panelSlider_Callback(hObject, ~)
        
        curSliderPosition = get(hObject, 'Value');
        set(subProcPanel, 'Position', subProcPanelPosition + [-unshownWidth*curSliderPosition/widthSubGUI, 0, 0, 0]);
        
    end
        
end

function generateVariablesPanelListing
    
    % Initialize Variables
    filetypeTextPosition = [panelBufferSpacing + panelBevelOffset + textBufferSpacing, heightGUI - panelBufferSpacing - panelBevelOffset - panelTitleHeights - textBufferSpacing - textIconHeight, widthGUI*experimentVariablesPanelWidthFraction - 3*textBufferSpacing - 2*panelBevelOffset - panelBufferSpacing - answerFieldDropDownWidth, textIconHeight];
    filetypeInputPosition = filetypeTextPosition + [filetypeTextPosition(3) + textBufferSpacing, 0, answerFieldDropDownWidth - filetypeTextPosition(3), 0];
    templateSizeTextPosition = filetypeTextPosition + [0, - textBufferSpacing - 2*textIconHeight, 0, textIconHeight];
    templateSizeInputPosition = templateSizeTextPosition + [templateSizeTextPosition(3) + textBufferSpacing + (filetypeInputPosition(3) - answerFieldEditWidth)/2, 0, answerFieldEditWidth - templateSizeTextPosition(3), inputIconHeight - textIconHeight];
    framerateTextPosition = templateSizeTextPosition + [0, - textBufferSpacing - 2*textIconHeight, 0, 0];
    framerateInputPosition = framerateTextPosition + [framerateTextPosition(3) + textBufferSpacing + (filetypeInputPosition(3) - answerFieldEditWidth)/2, 0, answerFieldEditWidth - framerateTextPosition(3), inputIconHeight - textIconHeight];
    scaleTextPosition = framerateTextPosition + [0, - textBufferSpacing - 2*textIconHeight, 0, 0];
    scaleInputPosition = scaleTextPosition + [scaleTextPosition(3) + textBufferSpacing + (filetypeInputPosition(3) - answerFieldEditWidth)/2, 0, answerFieldEditWidth - scaleTextPosition(3), inputIconHeight - textIconHeight];
    resReductionTextPosition = scaleTextPosition + [0, - textBufferSpacing - 2*textIconHeight, 0, 0];
    resReductionInputPosition = resReductionTextPosition + [resReductionTextPosition(3) + textBufferSpacing + (filetypeInputPosition(3) - answerFieldEditWidth)/2, 0, answerFieldEditWidth - resReductionTextPosition(3), inputIconHeight - textIconHeight];
    noVariableChangesPosition = resReductionTextPosition + [0, - textBufferSpacing - 3*textIconHeight, -resReductionTextPosition(3) + widthGUI*experimentVariablesPanelWidthFraction - 3*textBufferSpacing - 2*panelBevelOffset, textIconHeight];
    
    % Determine if user has modified anything and, if so, don't let them
    % modify the variables anymmore (analysis assumes original variables)
    shouldEnable = 'on';
    for i=1:size(currentAnalysisPerformed,2)
        for j=1:size(currentAnalysisPerformed(i).bools,1)
            if(sum(currentAnalysisPerformed(i).bools(j,1:4))>0)
                shouldEnable = 'off';
            end
        end
    end
    
    % Filetype text
    uicontrol('Parent',f,...
        'Style','text',...
        'String','Load which filetype?',...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',filetypeTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','popupmenu',...
        'String',{'*.tif', '*.png'},...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',filetypeInputPosition,...
        'Callback', {@filetype_Callback},...
        'Enable', 'off',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Template size text
    uicontrol('Parent',f,...
        'Style','text',...
        'String','Smallest template size for PIV tracking?',...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',templateSizeTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',analysisVariables{2},...
        'ForegroundColor',textFGColor,...
        'Position',templateSizeInputPosition,...
        'Callback', {@templateSize_Callback},...
        'Enable', shouldEnable,...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Framerate text
    uicontrol('Parent',f,...
        'Style','text',...
        'String','What is the framerate of the video (fps)?',...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',framerateTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',analysisVariables{3},...
        'ForegroundColor',textFGColor,...
        'Position',framerateInputPosition,...
        'Callback', {@framerate_Callback},...
        'Enable', shouldEnable,...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Scale text
    uicontrol('Parent',f,...
        'Style','text',...
        'String','What is the original spatial scale of the images (um/pix)?',...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',scaleTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',analysisVariables{4},...
        'ForegroundColor',textFGColor,...
        'Position',scaleInputPosition,...
        'Callback', {@scale_Callback},...
        'Enable', shouldEnable,...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Resolution reduction text, input
    uicontrol('Parent',f,...
        'Style','text',...
        'String','What (linear) factor would you like to reduce your image size by (1=none)?',...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',resReductionTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',analysisVariables{5},...
        'ForegroundColor',textFGColor,...
        'Position',resReductionInputPosition,...
        'Callback', {@resReduction_Callback},...
        'Enable', shouldEnable,...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % If user has run analysis, inform them they can no longer change
    % variables
    if(strcmp(shouldEnable,'off'))
        uicontrol('Parent',f,...
            'Style','text',...
            'String','Some Analysis is complete: Variables may no longer be changed',...
            'BackgroundColor',panelColor,...
            'ForegroundColor',[1,0,0],...
            'Position',noVariableChangesPosition,...
            'HorizontalAlignment','left',...
            'FontName','Gill Sans',...
            'FontSize',textSize);
    end
    
    function filetype_Callback(~, ~)
        
        msgbox('Currently not working');
        
    end
    
    function templateSize_Callback(hObject, ~)
        
        analysisVariables{2} = get(hObject,'String');
        
    end
    
    function framerate_Callback(hObject, ~)
        
        analysisVariables{3} = get(hObject,'String');
        
    end
    
    function scale_Callback(hObject, ~)
        
        analysisVariables{4} = get(hObject,'String');
        
    end
    
    function resReduction_Callback(hObject, ~)
        
        analysisVariables{5} = get(hObject,'String');
        
    end
    
end

function generateAnalysisPanelListing
    
    % Initialize Variables
    buttonWidth = 100;
    buttonHeight = 30;
    %textBufferSpacing = 4;
    fourierUpperBoundTextPosition = [analysisTitlePosition(1) + textBufferSpacing, analysisTitlePosition(2) - analysisTitlePosition(4) - panelBufferSpacing, widthGUI*experimentVariablesPanelWidthFraction - 3*textBufferSpacing - 2*panelBevelOffset - panelBufferSpacing - answerFieldEditWidth, 2*textIconHeight];
    fourierUpperBoundInputPosition = fourierUpperBoundTextPosition + [fourierUpperBoundTextPosition(3) + textBufferSpacing, 0, answerFieldEditWidth - fourierUpperBoundTextPosition(3), - textIconHeight];
    fourierLowerBoundTextPosition = fourierUpperBoundTextPosition + [0, - 2*textIconHeight - panelBufferSpacing, 0, 0];
    fourierLowerBoundInputPosition = fourierLowerBoundTextPosition + [fourierLowerBoundTextPosition(3) + textBufferSpacing, 0, answerFieldEditWidth - fourierLowerBoundTextPosition(3), - textIconHeight];
    VideoTimeTextPosition = fourierLowerBoundTextPosition + [0, - 2*textIconHeight - panelBufferSpacing, -fourierLowerBoundTextPosition(3)/4, 0];
    VideoTimeFirstInputPosition = VideoTimeTextPosition + [VideoTimeTextPosition(3) + textBufferSpacing, 0, answerFieldEditWidth/3 - VideoTimeTextPosition(3) + fourierLowerBoundTextPosition(3)/12, - textIconHeight];
    VideoXTextPosition = VideoTimeTextPosition + [0, - 2*textIconHeight - panelBufferSpacing, 0, 0];
    VideoXFirstInputPosition = VideoTimeFirstInputPosition + [0, - 2*textIconHeight - panelBufferSpacing, 0, 0];
    closeButtonPosition = [widthGUI*experimentVariablesPanelWidthFraction - panelBevelOffset - textBufferSpacing- buttonWidth, panelBufferSpacing + panelBevelOffset, buttonWidth, buttonHeight];
    analyzeButtonPosition = [panelBufferSpacing + panelBevelOffset + textBufferSpacing, panelBufferSpacing + panelBevelOffset, buttonWidth, buttonHeight];
    collectMotilityAnalysisPosition = [(analyzeButtonPosition(1) + closeButtonPosition(1))/2, (analyzeButtonPosition(2) + closeButtonPosition(2))/2, buttonWidth, buttonHeight];
    
    % Fourier upper bound text
    uicontrol('Parent',f,...
        'Style','text',...
        'String','Frequency upper bound on FFT (units min^-1)?',...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',fourierUpperBoundTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Fourier upper bound number
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',analysisVariables{6},...
        'ForegroundColor',textFGColor,...
        'Position',fourierUpperBoundInputPosition,...
        'Callback', {@fourierUpperBoundInput_Callback},...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Fourier lower bound text
    uicontrol('Parent',f,...
        'Style','text',...
        'String','Frequency lower bound on FFT (units min^-1)?',...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',fourierLowerBoundTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);

    % Fourier lower bound number
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',analysisVariables{7},...
        'ForegroundColor',textFGColor,...
        'Position',fourierLowerBoundInputPosition,...
        'Callback', {@fourierLowerBoundInput_Callback},...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Video StartT DeltaT EndT times percentages text
    uicontrol('Parent',f,...
        'Style','text',...
        'String',sprintf('Video Times:\nStart%% DeltaFrames End%%'),...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',VideoTimeTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Video StartT DeltaT EndT times percentages numbers
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',num2str(PIVVideoParams(1,1)),...
        'ForegroundColor',textFGColor,...
        'Position',VideoTimeFirstInputPosition,...
        'Callback', {@VideoParamsInput_Callback,1,1},...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',num2str(PIVVideoParams(1,2)),...
        'ForegroundColor',textFGColor,...
        'Position',VideoTimeFirstInputPosition + [answerFieldEditWidth/3 + fourierLowerBoundTextPosition(3)/12, 0, 0, 0],...
        'Callback', {@VideoParamsInput_Callback,1,2},...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',num2str(PIVVideoParams(1,3)),...
        'ForegroundColor',textFGColor,...
        'Position',VideoTimeFirstInputPosition + [2*(answerFieldEditWidth/3 + fourierLowerBoundTextPosition(3)/12), 0, 0, 0],...
        'Callback', {@VideoParamsInput_Callback,1,3},...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Video StartX DeltaX EndX X percentages text
    uicontrol('Parent',f,...
        'Style','text',...
        'String',sprintf('Video X Positions:\nStart%% ResReduction End%%'),...
        'BackgroundColor',panelColor,...
        'ForegroundColor',textFGColor,...
        'Position',VideoXTextPosition,...
        'HorizontalAlignment','left',...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Video StartX DeltaX EndX X percentages numbers
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',num2str(PIVVideoParams(2,1)),...
        'ForegroundColor',textFGColor,...
        'Position',VideoXFirstInputPosition,...
        'Callback', {@VideoParamsInput_Callback,2,1},...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',num2str(PIVVideoParams(2,2)),...
        'ForegroundColor',textFGColor,...
        'Position',VideoXFirstInputPosition + [answerFieldEditWidth/3 + fourierLowerBoundTextPosition(3)/12, 0, 0, 0],...
        'Callback', {@VideoParamsInput_Callback,2,2},...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    uicontrol('Parent',f,...
        'Style','edit',...
        'String',num2str(PIVVideoParams(2,3)),...
        'ForegroundColor',textFGColor,...
        'Position',VideoXFirstInputPosition + [2*(answerFieldEditWidth/3 + fourierLowerBoundTextPosition(3)/12), 0, 0, 0],...
        'Callback', {@VideoParamsInput_Callback,2,3},...
        'FontName','Gill Sans',...
        'FontSize',textSize);
    
    % Close button
    uicontrol('Parent',f,...
        'Style','pushbutton',...
        'String','Close Program',...
        'Position',closeButtonPosition,...
        'Callback',{@closeButton_Callback});
    
    % Collect motility button
    uicontrol('Parent',f,...
        'Style','pushbutton',...
        'String','Collect Analysis',...
        'Position',collectMotilityAnalysisPosition,...
        'Callback',{@collectMotilityAnalysis_Callback});
    
    % Analyze selection button
    uicontrol('Parent',f,...
        'Style','pushbutton',...
        'String','Analyze Selection',...
        'Position',analyzeButtonPosition,...
        'Callback',{@analyzeButton_Callback});
    
    % Fourier upper bound callback
    function fourierUpperBoundInput_Callback(hObject, ~)
        
        analysisVariables{6} = get(hObject,'String');
        
    end
    
    % Fourier lower bound callback
    function fourierLowerBoundInput_Callback(hObject, ~)
        
        analysisVariables{7} = get(hObject,'String');
        
    end
    
    % Functions for changing the PIVVideoParams
    function VideoParamsInput_Callback(hObject, ~, ii, jj)
        PIVVideoParams(ii,jj) = str2double(get(hObject,'String'));
    end
    
    % Close button callback
    function closeButton_Callback(~, ~)
        close all;
    end
    
    % Collect motility analysis button callback
    function collectMotilityAnalysis_Callback(~,~)
        
        motilityParams = collectMotilityAnalysis(mainAnalysisDirectory);
        assignin('base', 'motilityParams', motilityParams);
        
    end
    
    % Analyze motility button callback
    function analyzeButton_Callback(~, ~)
        % Close the program (it will open again with updated information)
        close all;
        
        % Perform PIV
        performPIV(mainExperimentDirectory, mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, rawPIVOutputName)
        
        % Obtain motility masks
        obtainMotilityMasks(mainExperimentDirectory, mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, maskFileOutputName)
        
        % Perform interpolation
        performMaskInterpolation(mainExperimentDirectory, mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, interpolationOutputName, rawPIVOutputName, maskFileOutputName)
        
        % Analyze Data
        performMotilityDataAnalysis(mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, motilityParametersOutputName, interpolationOutputName, GUISize)
        
        % Make PIV movie
        createAllChosenPIVMovies(mainExperimentDirectory, mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, PIVOutputName, PIVVideoParams, interpolationOutputName);
        
        % Reopen this program
        analyzeMotility(mainExperimentDirectory, mainAnalysisDirectory);
        
    end
    
end

end