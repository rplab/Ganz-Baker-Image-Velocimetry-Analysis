function fishParams = collectMotilityAnalysis( varargin )

%% Initialize variables

if( nargin < 1 )
        mainAnalysisDirectory = uigetdir(pwd, 'Main directory to contain/currently containing analysis');
    else
        mainAnalysisDirectory = varargin{ 1 };
end

[mainAnalysisDirectoryContents, mainAnalysisSubDirectoryContentsCell, nSubDirectories] = obtainDirectoryStructure(mainAnalysisDirectory);
currentAnalysisPerformedFile = load(strcat(mainAnalysisDirectory,filesep,'currentAnalysesPerformed.mat'));
currentAnalysisPerformed = currentAnalysisPerformedFile.currentAnalysisPerformed;
% disp('Warning: The function is not tested for accuracy since it was written on the fly. Consider writing your own');
index = 0;

%% Loop through all checked directories to generate PIV video
for i=1:nSubDirectories
    
    % Obtain the current directory size
    % nSubSubDirectories = size(mainAnalysisDirectoryContents(i).name, 1);
    subDire=dir(strcat(mainAnalysisDirectory,filesep,mainAnalysisDirectoryContents(i).name));
    subDire(strncmp({subDire.name}, '.', 1)) = []; % Removes . and ..
    subDire([subDire.isdir]==0) = []; % removes non-directories from list
    subFishDirect(i).name={subDire.name};
    nSFD=size(subFishDirect(i).name,2);
    curFolder = mainAnalysisDirectoryContents(i).name;
    % Loop through all checked subdirectories to perform PIV
    for j=1:nSFD
        
        if(currentAnalysisPerformed(i).bools(j,6))
            
            index = index + 1;
            currentAnalysisPerformed(i).directory
            % ObtainCurrentDirectory
            curAnDir = strcat(mainAnalysisDirectory, filesep, mainAnalysisDirectoryContents(i).name, filesep, mainAnalysisSubDirectoryContentsCell{1, i}(j).name);
            paramsFile = load(strcat(curAnDir, filesep, 'motilityParameters_Current.mat'));
            fishParams(index).Folder = curFolder;
            fishParams(index).SubFolder = mainAnalysisSubDirectoryContentsCell{1, i}(j).name;
            fishParams(index).FFTPowerPeak = paramsFile.fftPowerPeak;
            fishParams(index).FFTPeakFreq = paramsFile.fftPeakFreq;
            fishParams(index).FFTRPowerPeakSTD = paramsFile.fftRPowerPeakSTD;
            fishParams(index).FFTRPowerPeakMin = paramsFile.fftRPowerPeakMin;
            fishParams(index).FFTRPowerPeakMax = paramsFile.fftRPowerPeakMax;
            fishParams(index).WaveFrequency = paramsFile.waveFrequency;
            fishParams(index).WaveSpeedSlope = paramsFile.waveSpeedSlope;
            fishParams(index).BByFPS = paramsFile.BByFPS;
            fishParams(index).SigB = paramsFile.sigB;
            fishParams(index).WaveFitRSquared = paramsFile.waveFitRSquared;
            %        fishParams(index).XCorrMaxima = paramsFile.xCorrMaxima;
            fishParams(index).AnalyzedDeltaMarkersOne = paramsFile.analyzedDeltaMarkers(1);
            if(~isnan(paramsFile.analyzedDeltaMarkers(1)))
                fishParams(index).AnalyzedDeltaMarkersTwo = paramsFile.analyzedDeltaMarkers(2);
            else
                fishParams(index).AnalyzedDeltaMarkersTwo = NaN;
            end
            fishParams(index).WaveAverageWidth = paramsFile.waveAverageWidth;
            
        end

    end
    
end

temp_table = struct2table(fishParams);
writetable(temp_table,strcat(mainAnalysisDirectory, filesep, 'fishParams.csv'))

end

%% Auxiliary functions

% obtainDirectoryStructure returns directory structures given a directory
function [directoryContents, subDirectoryContentsCell, nSubDirectories] = obtainDirectoryStructure(directory)
    
    % Obtain main directory structure
    directoryContents = dir(directory); % Obtain all main directory contents
    directoryContents(~[directoryContents.isdir]) = []; % Remove non-directories
    directoryContents(strncmp({directoryContents.name}, '.', 1)) = []; % Removes . and .. and hidden files
    nSubDirectories = size(directoryContents, 1);
    subDirectoryContentsCell = cell(1, nSubDirectories);
    
    % Loop through all sub-directory contents to obtain contents
    if(nSubDirectories > 0)
        for i = 1:nSubDirectories
            
            % Obtain main directory structure
            subDirectoryContents = dir(strcat(directory, filesep, directoryContents(i).name)); % Obtain all sub-directory contents
            subDirectoryContents(~[subDirectoryContents.isdir]) = []; % Remove non-directories
            subDirectoryContents(strncmp({subDirectoryContents.name}, '.', 1)) = []; % Removes . and .. and hidden files
            subDirectoryContentsCell{i} = subDirectoryContents;
            
        end
    end
end