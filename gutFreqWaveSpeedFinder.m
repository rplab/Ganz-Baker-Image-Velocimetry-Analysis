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
% Modified May 13, 2023 (Raghuveer Parthasarathy; avoid filtering in time if
%           there aren't enough datapoints)
% Modified July 10-12, 2023 (Raghuveer Parthasarathy)
%          Allow use of median-x-subtracted cross-correlation image 
%          for identification of waves and determination
%          of frequency from the t=0 intercept. 
%          Avoid complicated (and troublesome) filtering, and instead fit a 1D
%          Gaussian to the XCorr(t) at each position to find the max. at
%          that position. Calls gaussfit1DMLE.m
% Last modified: July 12, 2023

function [waveFrequency, waveSpeedSlope, BByFPS, sigB, waveFitRSquared, ...
    xCorrMaxima, analyzedDeltaMarkers, XCorrFig] = ...
    gutFreqWaveSpeedFinder( gutMesh, trueXCorr, fps, scale )

% Find peristaltic frequency, wave speed from cross-correlation
NSeconds=90;
framesOfFirstNSeconds=1:NSeconds*fps;
numXCorrTimes = size(trueXCorr,1);
if NSeconds*fps > numXCorrTimes
    framesOfFirstNSeconds = 1:numXCorrTimes;
    fprintf('Duration of images < %.1f seconds\n', NSeconds); 
    fprintf('   Showing all frames of XCorr (%.1f seconds)\n',  numXCorrTimes/fps);
end

forYlim = min([NSeconds, numXCorrTimes/fps]);
% Plot (1) the cross-correlation image, and (2) the cross-correlation image 
% with each row's median subtracted to get rid of position-independent noise
% User can choose which to use:
xCorrIm1 = trueXCorr(framesOfFirstNSeconds,:);
xCorrIm2 = trueXCorr(framesOfFirstNSeconds,:) - median(trueXCorr(framesOfFirstNSeconds,:), 2);
g2 = figure; subplot(1,2,1)
plotCorrelationImage(xCorrIm1, forYlim);
subplot(1,2,2)
plotCorrelationImage(xCorrIm2, forYlim);
xCorrChoice=inputdlg({'Identify waves on the (1) original or (2) x-median-subtracted cross-corr. plot:'},...
    'XCorr',1,{'1'});
xCorrPlotChoice=round(str2double(xCorrChoice(1)));
close(g2)
XCorrFig = figure('name', 'Cross-correlation');
switch xCorrPlotChoice
    case 1
        plotCorrelationImage(xCorrIm1, forYlim);
        xCorrIm = xCorrIm1;
    case 2
        plotCorrelationImage(xCorrIm2, forYlim);
        xCorrIm = xCorrIm2;
    otherwise
        beep
        pause(1)
        disp('Bad choice; plotting and using original (not med. subtracted).')
        pause(1)
        xCorrIm = xCorrIm1;
        plotCorrelationImage(xCorrIm2, forYlim);
end

% Obtain rough line around maxima
roughFitEstimate = impoly( 'Closed', false );
rFEPoly = getPosition( roughFitEstimate );
deltaTimeUser=(rFEPoly(2,2)-rFEPoly(1,2));
deltaMarkerUser=rFEPoly(2,1)-rFEPoly(1,1);
translateMarkerNumToMicron=scale*round(mean(diff(squeeze(gutMesh(1,:,1,1)))));
slopeUser=deltaTimeUser/deltaMarkerUser; % Note: x is in units of marker numbers, not pixels or microns
interceptUser=rFEPoly(1,2)-slopeUser*rFEPoly(1,1);

% Ask for which ranges of x and t the user wants to search for maxima
last_x = min(size(xCorrIm, 2), 39);
xTDlgAns=inputdlg({'What range of time should be searched around (s)? +-', ...
    'What range of x should we use? First x is at marker number ', ...
    'What range of x should we use? Last x is at marker number '}, 'Title',1,{'8','1',num2str(last_x)});
timeAroundToSearchForMax=str2double(xTDlgAns(1));
markerNumStartX=str2double(xTDlgAns(2));
markerNumEndX=str2double(xTDlgAns(3));

% Find maxima: fit a Gaussian to XCorr(t) around the user-input line, at each x
xCorrMaxima = zeros(markerNumEndX - markerNumStartX + 1, 1);
for j=1:(markerNumEndX - markerNumStartX + 1)
    centerFrame = round(fps*(slopeUser*j+interceptUser));
    minFrame = round(centerFrame - timeAroundToSearchForMax*fps);
    maxFrame = round(centerFrame + timeAroundToSearchForMax*fps);
    xCorr_to_fit = xCorrIm(minFrame:maxFrame, markerNumStartX + j - 1);
    % Add offset to fit; we only care about the center
    [~, x0, ~, ~] = gaussfit1DMLE(xCorr_to_fit - min(xCorr_to_fit(:)) + 0.1);
    xCorrMaxima(j) = x0 + minFrame - 1; % frame at which XCorr is max. at this position j

end
xes = markerNumStartX:markerNumEndX;

figure(XCorrFig)
hold on;
plot(xes, (xCorrMaxima-1)/fps, 'x:', 'markersize', 6, 'color', 0.3*[1 1 1]);

[A, ~, B, sigB, waveFitRSquared, ~, ~] = fitline(xes, xCorrMaxima); 
plot(1, (A-1)/fps, 'ko', 'markersize', 12); % plot the intercept, i.e. the t=0 period

waveFrequency=60/((A-1)/fps); % Units of per minutes, -1 for indexing A at 1
waveSpeedSlope=fps*translateMarkerNumToMicron/B; % Units of um/sec ((frames/sec)*(micron/marker)/(frames/marker))
BByFPS = B/fps;

analyzedDeltaMarkers=[markerNumStartX, markerNumEndX];

end

function plotCorrelationImage(xCorrIm, forYlim)
% For showing the cross-correlation plot
    imshow(xCorrIm,[], 'InitialMagnification','fit', ...
        'YData', [0, forYlim], 'Xdata', [1, size(xCorrIm,2)], ...
        'border', 'loose');
    ylabel('Corr. Time (s)')
    xlabel('Corr X (marker)')
    set(gca,'YDir','normal')
    axis square;
    axis on;
    colormap('Jet');
end
