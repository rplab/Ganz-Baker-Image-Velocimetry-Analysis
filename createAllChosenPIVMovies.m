% Function which generates a movie containing the original image data
% superimposed with the velocity vectors. This video can optionally only be
% a subset of the full data.
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
%        - PIVOutputName: A string for naming the generated video. Should not
%            contain the directory.
%        - PIVVideoParams: A 2x3 array. A 1 in the first entry refers to
%            time to use parameters (start, delta, end) and a 2 refers to
%            the positions to use (start, delta, end).
%        - interpolationOutputName: String of the name of the file which 
%            contains the interpolated motility velocity vectors.

function createAllChosenPIVMovies(mainExperimentDirectory, mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, PIVOutputName, PIVVideoParams, interpolationOutputName)

%% Initialize variables
nDirectories = size(analysisToPerform, 2);
currentAnalysisFile = load(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName)); % WARNING: Do not change this variable name without changing the save string below
currentAnalysisPerformed = currentAnalysisFile.currentAnalysisPerformed; % WARNING: Don't change this variable name

% Progress bar
progtitle = sprintf('Preparing for video...');
progbar = waitbar(0, progtitle);  % will display progress

%% Loop through all checked directories to generate PIV video
for i=1:nDirectories
    
    % Progress bar update
    waitbar(i/nDirectories, progbar, ...
        sprintf('Generating PIV animation for folder %d of %d', i, nDirectories));
    
    % Obtain the current directory size
    nSubDirectories = size(analysisToPerform(i).bools, 1);
    
    % Loop through all checked subdirectories to perform PIV
    for j=1:nSubDirectories
        
        % If we want to analyze it, do so, else skip
        if(analysisToPerform(i).bools(j,5) && analysisToPerform(i).bools(j,6))
            
            % ObtainCurrentDirectory
            curAnDir = strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name);
            curExpDir = strcat(mainExperimentDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name);
            
            % Perform the PIV video generation on the current fish
            createPIVMovie(curAnDir, curExpDir, analysisVariables, PIVVideoParams, PIVOutputName, interpolationOutputName);
            
            % Update currentAnalysisPerformed
            currentAnalysisPerformed(i).bools(j,5) = true;
            
            % Save currentAnalysisPerformed.mat
            save(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName),'currentAnalysisPerformed','analysisVariables'); % WARNING: If currentAnalysisPerformed name is changed, you'll have to manually change this string IN MANY LOCATIONS!!!
            
        end
    end
    
end

close(progbar);

end