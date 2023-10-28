% createPIVMovie_standalone.m
% based on createPIVMovie.m (Ryan Baker)
%
% Function that superimposes velocity vectors from PIV onto source image frames,
% creating a PIV movie.
% Can run independently of the analyzeMotility GUI
% 
% Inputs
%     imagePath : folder containing images; if empty, dialog box
%     analysisParamsPath: path containing analysisParamsFile. 
%                    Leave blank for current directory. If analysisParamsFile
%                    is selected from a dialog box, ...Path will be set by
%                    this. Don't include final "\"
%     analysisParamsFile : file containing analysis results for this image set; 
%                    probably "currentAnalysesPerformed.mat"
%                    if empty, dialog box
%     PIVoutputFilePath : path containing PIVoutputFile
%                    Leave blank for current directory. If analysisParamsFile
%                    is selected from a dialog box, ...Path will be set by
%                    this. Don't include final "\"
%     PIVoutputFile : file containing PIV output, including gutMesh, 
%                     gutMeshVels, gutMeshVelsPCoords, thetas
%                     if empty, dialog box
%     PIVVideoParams : array of PIV parameters; default [0, 1, 100; 0, 1, 100];
%                    first 3 are start : deltaFrames : end Frame percent
%                    next 3 are positions.
%     PIVOutputName : output file name, default 'PIVAnimation.avi';
%     do_histEq: default false; if true, apply histogram equalization to the images.
%     velMultiple : multiply velocity by this for quiver length; formerly
%                    5; new default 20.
%
% Outputs
%    None.
%    Will write the PIV movie to PIVoutputFilePath, filename PIVOutputName
%
% Last modified: Oct. 10, 2023 -- Raghuveer Parthasarathy
%    - Avoid histogram equalization (optional)
%    - Change arrow color
%    - Change arrow size scaling factor
%    - The above are implemented as hard-coded parameters. (Not ideal, but
%      easy...)
%
% To do:
%    - Make a version of this that avoids the GUI
%
% Oct. 10, 2023 -- Raghuveer Parthasarathy
% Last modified: Oct. 10, 2023 -- Raghuveer Parthasarathy

function createPIVMovie_standalone(imagePath, analysisParamsPath, ...
    analysisParamsFile, PIVoutputFilePath, PIVoutputFile, PIVVideoParams, ...
    PIVOutputName, do_histEq, velMultiple)

if ~exist('imagePath', 'var') || isempty(imagePath)
    imagePath = uigetdir([], 'Folder containing images');
end
if ~exist('analysisParamsPath', 'var') || isempty(analysisParamsPath)
    analysisParamsPath = pwd;
end
if ~exist('analysisParamsFile', 'var') || isempty(analysisParamsFile)
    [analysisParamsFile, analysisParamsPath] = uigetfile('', 'Analysis results (currentAnalysesPerformed.mat)');
end
if ~exist('PIVoutputFilePath', 'var') || isempty(PIVoutputFilePath)
    PIVoutputFilePath = pwd;
end
if ~exist('PIVoutputFile', 'var') || isempty(PIVoutputFile)
    [PIVoutputFile, PIVoutputFilePath] = uigetfile('', 'PIV output file (processedPIVOutput_Current.mat)');
end
if ~exist('PIVVideoParams', 'var') || isempty(PIVVideoParams)
    PIVVideoParams = [0, 1, 100; 0, 1, 100];
end
if ~exist('PIVOutputName', 'var') || isempty(PIVOutputName)
    PIVOutputName = 'PIVAnimation.avi';
end
if ~exist('do_histEq', 'var') || isempty(do_histEq)
    do_histEq = false;
end
if ~exist('velMultiple', 'var') || isempty(velMultiple)
    velMultiple = 20;
end

load(strcat(analysisParamsPath,filesep,analysisParamsFile), 'analysisVariables');

quiverColor = [0.9 0.6 0.3]; % For Quiver

% % Load data and images from directory, define variables
suffix = analysisVariables{1};

%templateSize = str2double(analysisVariables{2});
%fps = str2double(analysisVariables{3});
%origMicronsPerPixel = str2double(analysisVariables{4});

origResReduction = str2double(analysisVariables{5});
% load(strcat(curAnDir,filesep,interpolationOutputName,'_Current.mat')); % Assumes file has gutMesh, gutMeshVels, gutMeshVelsPCoords, thetas
load(strcat(PIVoutputFilePath, filesep, PIVoutputFile))

startTP = PIVVideoParams(1,1)/100;
deltaF = PIVVideoParams(1,2);
endTP = PIVVideoParams(1,3)/100;
startX = PIVVideoParams(2,1)/100;
resReduceNew = PIVVideoParams(2,2);
endX = PIVVideoParams(2,3)/100;
resReduce = resReduceNew*origResReduction;

direc = dir([imagePath,filesep,suffix]); 
baseFilenames={};
[baseFilenames{1:length(direc),1}] = deal(direc.name);
baseFilenames = sortrows(baseFilenames); %sort all image files
amount = length(baseFilenames);
count=1; % Linear index that travels through all multipages monotonically
filenames = {};

writerObj = VideoWriter(strcat(PIVoutputFilePath,filesep,PIVOutputName),'Uncompressed AVI');
showOnlyHorizontalComponent = true;

% Progress bar
progtitle1 = sprintf('Obtaining image information...');
progbar1 = waitbar(0, progtitle1);  % will display progress

% Build a single large array of filenames
if(strcmp(suffix, '*.tif'))
    
    % If we have tifs, they may be multipage, so build a new name that
    % includes the appropriate index
    for i=1:amount
        
        % Progress bar update
        waitbar(i/amount, progbar1, ...
            sprintf('Obtaining image information for image stack %d of %d', i, amount));
        
        info = imfinfo([imagePath filesep baseFilenames{i}]);
        nI=size(info,1);
        
        for j=1:nI
            
            filenames{count}.name = info.Filename; %#ok since I can't think of an easy/smart way of preallocating
            filenames{count}.index = j; %#ok since I can't think of an easy/smart way of preallocating
            
            count = count+1;
            
        end
        
    end
    
else
       
    % If we don't have tifs, then the filenames is the array we want
    for i=1:amount
        
        filenames{i}.name = baseFilenames{i}; %#ok since I can't think of an easy/smart way of preallocating
        
    end
    
end

close(progbar1);

% Video properties
nF=size(filenames,2);
info = imfinfo([imagePath filesep baseFilenames{1}]);
numRows = info(1).Height;
numCols = info(1).Width;

% Reduced video properties
startingFrame = round(nF*startTP + 1);
endingFrame = round(nF*endTP);
startingPosition = round(numCols*startX + 1);
endingPosition = round(numCols*endX);

if do_histEq
    firstImHistExample = histeq(imread(fullfile(filenames{1}.name), 'Index', filenames{1}.index,'PixelRegion', {[1 resReduce numRows], [startingPosition resReduce endingPosition]})); % read images
    maxI = double(max(firstImHistExample(:)));
    minI = double(min(firstImHistExample(:)));
end

% Open video writing code, initialize settings
open(writerObj);
figure;
set(gcf,'Renderer','zbuffer');
% totalFrames = ceil((endingFrame - startingFrame + 1)/deltaF);

for i=startingFrame:deltaF:endingFrame - 1
    
    % Get image
    %imageIndex = (i - startingFrame + deltaF)/deltaF - (i~=1); % The logical subtraction is because gutMeshVels and others only have N - 1 frames from differencing; just make the first frame the same as the second
    curImage = imread(fullfile(filenames{i}.name), 'Index', filenames{i}.index, 'PixelRegion', {[1 resReduce numRows], [startingPosition resReduce endingPosition]});
    if do_histEq
        im = (double(imhistmatch(curImage,firstImHistExample)) - minI)/(maxI - minI); % read images, histEQ them
    else
        % no equalization
        im = double(curImage);
    end
   
    % Get full vector field for quiver plot
    curIndex = i - (i~=1); % The logical subtraction is because gutMeshVels and others only have N - 1 frames from differencing; just make the first frame the same as the second
    qx=gutMesh(:,:,1);
    qy=gutMesh(:,:,2);
    qu=velMultiple*gutMeshVels(:,:,1,curIndex);
    qv=velMultiple*gutMeshVels(:,:,2,curIndex);
    
    % Get local representations of components
    qup=velMultiple*squeeze(gutMeshVelsPCoords(:,:,1,curIndex));
    qvp=velMultiple*squeeze(gutMeshVelsPCoords(:,:,2,curIndex));
    % Re-project each component onto x y axes for visual representation
    qupx=qup; % This is just to get the size right
    qupx(1:end/2,:)=qup(1:end/2,:)*cos(thetas(1)); % Top
    qupx((end/2+1):end,:)=qup((end/2+1):end,:)*cos(thetas(2)); % Bottom
    qupy=qup;
    qupy(1:end/2,:)=qup(1:end/2,:)*sin(thetas(1));
    qupy((end/2+1):end,:)=qup((end/2+1):end,:)*sin(thetas(2));
    qvpx=qvp;
    qvpx(1:end/2,:)=-qvp(1:end/2,:)*sin(thetas(1));
    qvpx((end/2+1):end,:)=-qvp((end/2+1):end,:)*sin(thetas(2));
    qvpy=qvp;
    qvpy(1:end/2,:)=qvp(1:end/2,:)*cos(thetas(1));
    qvpy((end/2+1):end,:)=qvp((end/2+1):end,:)*cos(thetas(2));
    
    if(showOnlyHorizontalComponent)
        qv = qv.*0;
    end
    
    imshow(im,[]);
    hold on;
    quiver(qx,qy,qu,qv,0,'Color', quiverColor,'LineWidth',2);
    %quiver(qx,qy,qupx,qupy,0,'b');
    %quiver(qx,qy,qvpx,qvpy,0,'b');
    hold off;
    
    % Write image to video
    frame = getframe;
    writeVideo(writerObj,frame);
    
end

% Close video writing code
close(writerObj);

% % Load data and images from directory
% matData=dir(strcat(imPath,filesep,suffix)); % Assumes filesep at end of imPath
% load(strcat(imPath,filesep,matData(1).name));
% ims=dir(strcat(imPath,filesep,filetype));
% nT=size(ims,1)-1; % differencing frames obviously leads to n-1 frames

end