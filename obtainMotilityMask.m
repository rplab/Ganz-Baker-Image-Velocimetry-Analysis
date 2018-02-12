% Function which shows a histogram equalized version of the first image in
% a subdirectory and allows the user to draw a mask over it by selecting
% vertices and drawing a polygon. Double clicking closes the polygon, with
% the inside being the part which is analyzed. Then the user draws a center
% line. Double clicking ends the centerline drawing.
%
% Inputs:- imPath: A string representing the path to the saved images in
%            question.
%        - filetype: A string (with wildcard) representation of the
%            filetype of the images.
%
% Outputs:- gutOutline: A handle to the polygon object. This output is
%             no longer used.
%         - gutOutlinePoly: A 2xN representation of all the vertices (in image
%             pixel coordinates) used to represent the mask. The last entry
%             and the first entry are assumed to be connected.
%         - gutMiddleTop: A handle to the polygon object for the center
%             line. This output is no longer used.
%         - gutMiddleBottom: Same as gutMiddleTop. Used to contain the
%             location of the bottom edge of the luman.
%         - gutMiddlePolyTop: A 2xM representation of all the vertices (in
%             image pixel coordinates) used to represent the center line.
%             The last entry is not connected to the first entry.
%         - gutMiddlePolyBottom: Same as gutMiddlePolyTop.


function [gutOutline, gutOutlinePoly, gutMiddleTop, gutMiddleBottom, gutMiddlePolyTop, gutMiddlePolyBottom] = obtainMotilityMask(imPath,filetype,resReduce)

%% Read first image
ims=dir(strcat(imPath,filesep,filetype));
info = imfinfo(strcat(imPath,filesep,ims(1).name));
nNumCols = info(1).Height;
nNumRows = info(1).Width;
im=imread(strcat(imPath,filesep,ims(1).name),'PixelRegion', {[1 resReduce nNumCols], [1 resReduce nNumRows]});

% Initialize variables
continueBool = 0;

%% Draw gut outline, gut center
while (continueBool ~= 1)
    
    close all;
    
    % Draw the boundary and gut center. Try to draw the centerline in a way
    % you'd want orthogonal vectors to orient themselves
    figure;
    imH = imshow(histeq(im), []);
    %imC = imcontrast(imH);
    gutOutline = impoly;
    gutOutlinePoly = getPosition(gutOutline);
    setColor(gutOutline, 'red'); 
    gutMiddleTop = impoly('Closed', false);
    gutMiddleBottom = gutMiddleTop;
    gutMiddlePolyTop=getPosition(gutMiddleTop);
    gutMiddlePolyBottom = getPosition(gutMiddleBottom);
    setColor(gutMiddleTop,'green'); 
    setColor(gutMiddleBottom,'green'); 
    
    % Prompt user to see if they want to continue or redraw
    continueBool=menu('Are you satisfied with your drawing?', ...
        'Yes', 'No');
    
    close all;

end

close all;

end