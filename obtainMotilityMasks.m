% Function which loops through each directory and subdirectory and, if the
% user wants to create a mask and continue the analysis, prompts them with
% a image onto which they can draw vertices representative of that mask.
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
%        - maskFileOutputName: String of the name of the file which will 
%            contain the mask vertices.

function obtainMotilityMasks(mainExperimentDirectory, mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, maskFileOutputName)

%% Initialize variables
nDirectories = size(analysisToPerform, 2);
currentAnalysisFile = load(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName)); % WARNING: Do not change this variable name without changing the save string below
currentAnalysisPerformed = currentAnalysisFile.currentAnalysisPerformed; % WARNING: Don't change this variable name

%% Loop through all checked directories to obtain image masks
for i=1:nDirectories
    
    % Obtain the current directory size
    nSubDirectories = size(analysisToPerform(i).bools, 1);
    
    % Loop through all checked subdirectories to perform PIV
    for j=1:nSubDirectories
        
        % If we want to analyze it, do so, else skip
        if(analysisToPerform(i).bools(j,2) && analysisToPerform(i).bools(j,6))
            
            % ObtainCurrentDirectory
            curDir = strcat(mainExperimentDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name);
            
            % Perform mask creation
            [gutOutline, gutOutlinePoly, gutMiddleTop, gutMiddleBottom, gutMiddlePolyTop, gutMiddlePolyBottom] = obtainMotilityMask(curDir, analysisVariables{1}, str2double(analysisVariables{5})); %#ok since it is saved WARNING: Don't change these variable names
            
            % Save <maskFileOutputName>_Current.mat, <maskFileOutputName>_<date>.mat
            save(strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name, filesep, maskFileOutputName, '_Current'), 'gutOutline', 'gutOutlinePoly', 'gutMiddleTop', 'gutMiddleBottom', 'gutMiddlePolyTop', 'gutMiddlePolyBottom');
            save(strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name, filesep, maskFileOutputName, '_', date), 'gutOutline', 'gutOutlinePoly', 'gutMiddleTop', 'gutMiddleBottom', 'gutMiddlePolyTop', 'gutMiddlePolyBottom');
            
            % Update currentAnalysisPerformed
            currentAnalysisPerformed(i).bools(j,2) = true;
            
            % Save currentAnalysisPerformed.mat
            save(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName),'currentAnalysisPerformed','analysisVariables'); % WARNING: If currentAnalysisPerformed name is changed, you'll have to manually change this string IN MANY LOCATIONS!!!
            
        end
    end
    
end

end