% Function which shows an image of the cross correlation, prompts the user
% for a line showing the program roughly where first set of correlation
% maxima are, and then processes that subset of the data (filtering) 
% to extract wave speed, frequency, variation, etc.
%
% Inputs:- gutMesh: A NxM array of vertex locations for a grid morphed to 
%             fit inside of the mask. The slopes for each column of 
%             vertices is given by thetas.
%        - trueXCorr: The cross-correlation surface generated elsewhere,
%             with flips taken into account.
%        - fps: Frames per second of the movie being analyzed. 
%        - scale: The reduced scale of the images
%          (resolution*resolutionReduction), in units of microns/pixel.
%
% Outputs:- waveFrequency: A number representing what we now call Frequency.
%             Obtained from the XCorr. In units of per minutes.
%         - waveSpeedSlope: A number representing the slope of the
%             cross-correlation peaks, and thus the speed of a wave as it
%             travels down the gut. In units of microns/second.
%         - BByFPS: Inverse wave speed slope. In units of seconds per
%             pixel.
%         - sigB: Slope uncertainty for the line fit to the
%             cross-correlation maxima.
%         - waveFitRSquared: R-squared value for the line fit to the
%             cross-correlation maxima.
%         - xCorrMaxima: A vector containing the time at which the
%             cross-correlation is maximum for each delta x.
%         - analyzedDeltaMarkers: a 1x2 vector containing the start and
%             finish of the delta x's to use when fitting a line to the
%             cross-correlation maxima.
%         - g: A handle to the figure that is generated. Useful for
%             shutting down by command later.

% To do:
%    -- Make "NSeconds" not be hard-coded to 90. Use all the timepoints, 
%       or have this be a more transparent parameter somewhere. [RP Sept 17, 2022]
%    -- .
%
% Ryan P. Baker
% Modified: Sept. 9, 2018 (Raghuveer Parthasarathy; minor change to display of Xcorr)
% Modified: Sept. 17, 2022 (Raghuveer Parthasarathy; rounding of dTime, line 90)
%           Modify filter settings to not be so crudely hardcoded.
%           Check that filter creation doesn't give an error; expand if needed.
% Last modified: Sept. 17, 2022

function [waveFrequency, waveSpeedSlope, BByFPS, sigB, waveFitRSquared, ...
    xCorrMaxima, analyzedDeltaMarkers, g] = ...
    gutFreqWaveSpeedFinder( gutMesh, trueXCorr, fps, scale )

% Find peristaltic frequency, wave speed from cross-correlation
g = figure('name', 'Cross-correlation');
NSeconds=90;
framesOfFirstNSeconds=1:NSeconds*fps;
numXCorrTimes = size(trueXCorr,1);
if NSeconds*fps > numXCorrTimes
    framesOfFirstNSeconds = 1:numXCorrTimes;
    fprintf('Duration of images < %.1f seconds\n', NSeconds); 
    fprintf('   Showing all frames of XCorr (%.1f seconds)\n',  numXCorrTimes/fps);
end
imshow(trueXCorr(framesOfFirstNSeconds,:),[], 'InitialMagnification','fit', ...
    'YData', [0, min([NSeconds, numXCorrTimes/fps])], 'Xdata', [1, size(trueXCorr,2)], ...
    'border', 'loose');
ylabel('Corr. Time (s)')
xlabel('Corr X (marker)')
set(gca,'YDir','normal')
axis square;
axis on;
colormap('Jet');

% Obtain rough line around maxima
roughFitEstimate = impoly( 'Closed', false );
rFEPoly = getPosition( roughFitEstimate );
deltaTimeUser=(rFEPoly(2,2)-rFEPoly(1,2));
deltaMarkerUser=rFEPoly(2,1)-rFEPoly(1,1);
translateMarkerNumToMicron=scale*round(mean(diff(squeeze(gutMesh(1,:,1,1)))));
slopeUser=deltaTimeUser/deltaMarkerUser; % Note: x is in units of marker numbers, not pixels or microns
interceptUser=rFEPoly(1,2)-slopeUser*rFEPoly(1,1);

% Ask for which ranges of x and t the user wants to search rather than
% using initialized variables, smoothing size
last_x = min(size(trueXCorr(framesOfFirstNSeconds,:), 2), 39);
xTDlgAns=inputdlg({'What range of time should be searched around (s)? +-', ...
    'What range of x should we use? First x is at marker number ', ...
    'What range of x should we use? Last x is at marker number '}, 'Title',1,{'5','1',num2str(last_x)});
timeAroundToSearchForMax=str2double(xTDlgAns(1));
markerNumStartFreq=str2double(xTDlgAns(2));
markerNumEndFreq=str2double(xTDlgAns(3));

% Perform zero-phase digital filter to data
dMarker=(markerNumEndFreq-markerNumStartFreq);
%dTime=2*floor(fps*timeAroundToSearchForMax); % Always even
dTimeMin=round(min(slopeUser*(markerNumStartFreq:markerNumEndFreq)+interceptUser-timeAroundToSearchForMax));
dTimeMax=round(max(slopeUser*(markerNumStartFreq:markerNumEndFreq)+interceptUser+timeAroundToSearchForMax));
% Check that dTimeMin is positive (and at least 1/fps)

% This next if statement is needed for my filter, minimum size of 27 (using literals here is not wise, I know)
dTime=floor(dTimeMax-dTimeMin);
if(dTime*fps<27)
    addToDTimes = ceil((dTime*fps-27)/2);
    dTimeMin = dTimeMin - addToDTimes;
    dTimeMax = dTimeMax + addToDTimes;
    disp('Warning: Time range extended to allow filter to work. See gutFreqWaveSpeedFinder line 95 for more info');
end

if dTimeMin < 1/fps
    fprintf('gutFreqWaveSpeedFinder: dTimeMin is negative or too small (%.1f). Force 1/fps.\n', dTimeMin)
    dTimeMinOffset = 1/fps - dTimeMin;
    dTimeMin = 1/fps;
    dTimeMax = dTimeMax + dTimeMinOffset;
    fprintf('Also offset dTimeMax: now %0.1f.\n', dTimeMax)
end
dTime=floor(dTimeMax-dTimeMin);

% Make 'PassbandFrequency',  and 'StopbandFrequency' depend on fps. Note
% that values of 0.15 and 0.65 seem fine for 5 fps, so scale by this; max 1
% (normalized frequency)
LPFilt=designfilt('lowpassfir', 'PassbandFrequency', min(0.15*5/fps, 1.0), ...
        'StopbandFrequency', min(0.65*5/fps, 1.0), 'PassbandRipple', 1, ...
        'StopbandAttenuation', 60);
nCoefficients = length(LPFilt.Coefficients); % number of coefficients in the filter

% If the number of timepoints is not > 3*nCoefficents, filtering will not
% work. Adjust dTimeMax
if (dTime*fps + 1) < 3*nCoefficients
    fprintf('gutFreqWaveSpeedFinder: Insufficient timepoints. dTime*fps = %.1f; nCoefficients = %d\n', ...
        dTime*fps, nCoefficients)
    dTimeMax = (3*nCoefficients+1)/fps + dTimeMin;
    if dTimeMax*fps > size(trueXCorr,1)
        errordlg('gutFreqWaveSpeedFinder: insufficient timepoints for filter; can''t extend.')
    end
    dTime=floor(dTimeMax-dTimeMin);
    fprintf('   Extending dTimeMax; now: %.1f\n', dTimeMax);
end

reducedSmoothedVelocityMap=zeros(dMarker,dTime+1);
for i=markerNumStartFreq:markerNumEndFreq
    %tAroundToSearch=round(fps*(slopeUser*i+interceptUser));
    subsetTrueXCorr=squeeze(trueXCorr(fps*dTimeMin:fps*dTimeMax,i));
    reducedSmoothedVelocityMap(i,:)=filtfilt(LPFilt,subsetTrueXCorr);
end

% Find maxima
[~, xCorrMaxima]=max(reducedSmoothedVelocityMap,[],2); % Name is misleading, should be xCorrMaximaTimes I think?
xes = markerNumStartFreq:markerNumEndFreq;
[A, ~, B, sigB, waveFitRSquared, ~, ~] = fitline( xes, xCorrMaxima' );
waveFrequency=60/((A-1)/fps+dTimeMin); % Units of per minutes, -1 for indexing A at 1
waveSpeedSlope=fps*translateMarkerNumToMicron/B; % Units of um/sec ((frames/sec)*(micron/marker)/(frames/marker))
BByFPS = B/fps;

% % Fit to line, get slope/intercept for wave speed/frequency, variance about linear fit
% linearCoefs=polyfit(1:size(xCorrMaxima,1),xCorrMaxima',1);
% %fprintf('Intercept: %d \n',(linearCoefs(2)+dTimeMin-1));
% 
% % Find R^2
% yfit=polyval(linearCoefs,1:size(xCorrMaxima,1));
% yresid=xCorrMaxima'-yfit;
% SSresid = sum(yresid.^2); % Awful units!
% SStotal = (size(xCorrMaxima,1)-1) * var(xCorrMaxima);
% waveFitRSquared = 1 - SSresid/SStotal;

analyzedDeltaMarkers=[markerNumStartFreq, markerNumEndFreq];

end