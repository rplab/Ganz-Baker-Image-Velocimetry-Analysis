% Function which loops through each directory and subdirectory and, if the
% user wants to perform the analysis and use the data, uses the command
% line version of PIVLab to obtain the velocity vector field representative
% of the underlying image data sets.
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
%        - rawPIVOutputName: String of the name of the file which will 
%            contain the raw PIV vector information

function performPIV(mainExperimentDirectory, mainExperimentDirectoryContents, mainExperimentSubDirectoryContentsCell, mainAnalysisDirectory, analysisToPerform, analysisVariables, currentAnalysesPerformedFileName, rawPIVOutputName)

%% Initialize variables
nDirectories = size(analysisToPerform, 2);
currentAnalysisFile = load(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName)); % WARNING: Do not change this variable name without changing the save string below
currentAnalysisPerformed = currentAnalysisFile.currentAnalysisPerformed; % WARNING: Don't change this variable name

% Progress bar
progtitle = sprintf('Initializing PIV analysis...');
progbar = waitbar(0, progtitle);  % will display progress

%% Loop through all checked directories to perform PIV
for i=1:nDirectories
    
    % Obtain the current directory size
    nSubDirectories = size(analysisToPerform(i).bools, 1);
    
    % Progress bar update
    waitbar(i/nDirectories, progbar, ...
        sprintf('Obtaining PIV vectors for folder %d of %d', i, nDirectories));
    
    % Loop through all checked subdirectories to perform PIV
    for j=1:nSubDirectories
        
        % If we want to analyze it, do so, else skip
        if(analysisToPerform(i).bools(j,1) && analysisToPerform(i).bools(j,6))
            
            % ObtainCurrentDirectory
            curDir = strcat(mainExperimentDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name);
            
            % Perform PIV, obtain rawPIVOutput
            [p, s, x, y, u, v, typevector, imageDirectory, filenames, u_filt, v_filt, typevector_filt] = obtainRawPIVOutput(curDir, analysisVariables); %#ok since it is saved WARNING: Don't change these variable names
            
            % Save <rawPIVOutputName>_Current.mat, <rawPIVOutputName>_<date>.mat
            save(strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name, filesep, rawPIVOutputName, '_Current'), 'p', 's', 'x', 'y', 'u', 'v', 'typevector', 'imageDirectory', 'filenames', 'u_filt', 'v_filt', 'typevector_filt');
            save(strcat(mainAnalysisDirectory, filesep, mainExperimentDirectoryContents(i).name, filesep, mainExperimentSubDirectoryContentsCell{1, i}(j).name, filesep, rawPIVOutputName, '_', date), 'p', 's', 'x', 'y', 'u', 'v', 'typevector', 'imageDirectory', 'filenames', 'u_filt', 'v_filt', 'typevector_filt');
            
            % Update currentAnalysisPerformed
            currentAnalysisPerformed(i).bools(j,1) = true;
            
            % Save currentAnalysisPerformed.mat
            save(strcat(mainAnalysisDirectory, filesep, currentAnalysesPerformedFileName),'currentAnalysisPerformed','analysisVariables'); % WARNING: If currentAnalysisPerformed name is changed, you'll have to manually change this string IN MANY LOCATIONS!!!
            
        end
    end
    
end

close(progbar);

end