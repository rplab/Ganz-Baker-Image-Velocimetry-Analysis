% Function which loops through each directory and subdirectory and, if the
% user wants to interpolate the PIV output to fit inside their mask and
% continue the analysis, does exactly that.
%
% Inputs:- mainExperimentDirectory: Directory containing the raw
%            image data. Input 0 if you don't want to use it. Prompts for
%            directory if one isn't given. Push cancel if you don't want to
%            use it.
%        - mainExperimentDirectoryContents: A list of all main folder names
%            and other metadata.
%        - mainExperimentSubDirectoryContentsCell: A cell array in which
%            each element, representative of one main folder (see above),
%            contains a list of all subdirectories' names and other
%            metadata in that main folder.
%        - mainAnalysisDirectory: Directory that will/does contain the 
%            analyzed data. Prompts for directory if one isn't given.
%        - analysisToPerform: A structure holding which analysis should be
%            performed as a bool array.
%        - analysisVariables: A cell array of numbers and strings
%            containing the user variables specified in the GUI (top left).
%        - currentAnalysesPerformedFileName: String of the name of the file
%            which will contain a memory of which analysis has been
%            performed.
%        - interpolationOutputName: String of the name of the file which 
%            will contain the interpolated motility velocity vectors.
%        - rawPIVOutputName: String of the name of the file which contains
%            the raw PIV vector information to load.
%        - maskFileOutputName: String of the name of the file which 
%            contains the mask vertices to load.

function performMaskInterpolation(mainExperimentDirectory, mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, interpolationOutputName, rawPIVOutputName, maskFileOutputName)

%% Initialize variables
nDirectories = size(analysisToPerform, 2);
currentAnalysisFile = load(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName)); % WARNING: Do not change this variable name without changing the save string below
currentAnalysisPerformed = currentAnalysisFile.currentAnalysisPerformed; % WARNING: Don't change this variable name

% Progress bar
progtitle = sprintf('Preparing for interpolation...');
progbar = waitbar(0, progtitle);  % will display progress

%% Loop through all checked directories to perform mask interpolation
for i=1:nDirectories
    
    % Progress bar update
    waitbar(i/nDirectories, progbar, ...
        sprintf('Performing interpolation for folder %d of %d', i, nDirectories));
    
    % Obtain the current directory size
    nSubDirectories = size(analysisToPerform(i).bools, 1);
    
    % Loop through all checked subdirectories to perform interpolation
    for j=1:nSubDirectories
        
        % If we want to analyze it, do so, else skip
        if(analysisToPerform(i).bools(j,3) && analysisToPerform(i).bools(j,6))
            
            % ObtainCurrentDirectory
            curDir = strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name);
            expDir = strcat(mainExperimentDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name);
            
            % Perform mask creation
            [gutMesh, mSlopes, gutMeshVels, gutMeshVelsPCoords, thetas] = interpolatePIVVectorsInMask(curDir, expDir, analysisVariables{1}, str2double(analysisVariables{5}), rawPIVOutputName, maskFileOutputName); %#ok since it is saved WARNING: Don't change these variable names
            
            % Save <interpolationOutputName>_Current.mat, <interpolationOutputName>_<date>.mat
            save(strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name, filesep, interpolationOutputName, '_Current'), 'gutMesh','mSlopes','gutMeshVels','gutMeshVelsPCoords','thetas');
            save(strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name, filesep, interpolationOutputName, '_', date), 'gutMesh','mSlopes','gutMeshVels','gutMeshVelsPCoords','thetas');
            
            % Update currentAnalysisPerformed
            currentAnalysisPerformed(i).bools(j,3) = true;
            
            % Save currentAnalysisPerformed.mat
            save(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName),'currentAnalysisPerformed','analysisVariables'); % WARNING: If currentAnalysisPerformed name is changed, you'll have to manually change this string IN MANY LOCATIONS!!!
            
        end
    end
    
end

close(progbar);

end