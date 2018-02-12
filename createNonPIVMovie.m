% Function which...
%
% Remember to do histogram normalization

function createNonPIVMovie(imageDirectory, PIVVideoParams)

%% Create list of images inside specified directory
suffix = '*.tif';
startTP = PIVVideoParams(1,1)/100;
deltaF = PIVVideoParams(1,2);
endTP = PIVVideoParams(1,3)/100;
startX = PIVVideoParams(2,1)/100;
resReduce = PIVVideoParams(2,2);
endX = PIVVideoParams(2,3)/100;
direc = dir([imageDirectory,filesep,suffix]); 
baseFilenames={};
[baseFilenames{1:length(direc),1}] = deal(direc.name);
baseFilenames = sortrows(baseFilenames); %sort all image files
amount = length(baseFilenames);
count=1; % Linear index that travels through all multipages monotonically
filenames = {};

% Progress bar
progtitle1 = sprintf('Obtaining image data...');
progbar1 = waitbar(0, progtitle1);  % will display progress

% Build a single large array of filenames
if(strcmp(suffix, '*.tif'))
    
    % If we have tifs, they may be multipage, so build a new name that
    % includes the appropriate index
    for i=1:amount
        
        % Progress bar update
        waitbar(i/amount, progbar1, ...
            sprintf('Obtaining image data for image stack %d of %d', i, amount));
        
        info = imfinfo([imageDirectory filesep baseFilenames{i}]);
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
info = imfinfo([imageDirectory filesep baseFilenames{1}]);
numRows = info(1).Height;
numCols = info(1).Width;

% Reduced video properties
startingFrame = round(nF*startTP + 1);
endingFrame = round(nF*endTP);
startingPosition = round(numCols*startX + 1);
endingPosition = round(numCols*endX);

% Initialize video data (CAREFUL: This can be huge if not done correctly)
imStack = zeros(ceil(numRows/resReduce),ceil((endingPosition - startingPosition + 1)/resReduce), ceil((endingFrame - startingFrame + 1)/deltaF)); % You'd normally use the floor function rather than ceil, but this would lead to a classic "fencepost" error; we need to include the ends, not just the delta
firstImHistExample = histeq(imread(fullfile(filenames{1}.name), 'Index', filenames{1}.index,'PixelRegion', {[1 resReduce numRows], [startingPosition resReduce endingPosition]})); % read images
maxI = double(max(firstImHistExample(:)));
minI = double(min(firstImHistExample(:)));

%% Load video data into program

% Progress bar
progtitle = sprintf('Preparing to load video...');
progbar = waitbar(0, progtitle);  % will display progress

for i=startingFrame:deltaF:endingFrame
    
    imageIndex = (i - startingFrame + deltaF)/deltaF;
    totalFrames = ceil((endingFrame - startingFrame + 1)/deltaF);
    
    % Progress bar update
    waitbar(imageIndex/totalFrames, progbar, ...
        sprintf('Loading frame %d of %d', imageIndex, totalFrames));
    
    curImage = imread(fullfile(filenames{i}.name), 'Index', filenames{i}.index,'PixelRegion', {[1 resReduce numRows], [startingPosition resReduce endingPosition]});
    imStack(:,:,imageIndex) = (double(imhistmatch(curImage,firstImHistExample)) - minI)/(maxI - minI); % read images, histEQ them

end

close(progbar);

%% Show video
implay(imStack);

end