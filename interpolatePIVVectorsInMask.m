% Function which takes the raw PIV vectors (previously obtained) and
% interpolates them into a NxM grid morphed into the mask (previously
% obtained). It also projects these velocities along a new, non-euclidean
% coordinate system for longitudinal and transverse projections.
%
% Inputs:- curDir: The directory to obtain the raw PIV and mask data from.
%        - expDir: The directory containing the image data.
%        - imageType: A string of the filetype to load (including wildcard)
%        - resReduce: A integer representing what factor to reduce the
%            resolution by (e.g. 2).
%        - rawPIVOutputName: String of the name of the file which contains
%            the raw PIV vector information to load.
%        - maskFileOutputName: String of the name of the file which 
%            contains the mask vertices to load.
%
% Outputs:- gutMesh: A NxM array of vertex locations for a grid morphed to 
%             fit inside of the mask. The slopes for each column of 
%             vertices is given by thetas.
%         - mSlopes: A Mx2x2 vector of differences for each column. Each dx
%             and dy is calculated as the perpendicular to the center line. 
%             The first input is which column, the second is either dx or
%             dy, and the third input is top or bottom (which are now the
%             same thing).
%         - gutMeshVels: A NxMx2xT array of doubles representing
%             velocities. The first and second input correspond to which
%             vertex (see gutMesh), the third input corresponds to which
%             component of the velocity (x or y) and the last input
%             corresponds to which frame.
%         - gutMeshVelsPCoords: A NxMx2xT array of doubles representing
%             velocities. These velocities are in a new primed coordinate 
%             system, longitudinal or transverse. Note that the space, as a
%             result, is non-euclidean. The first and second input 
%             correspond to which vertex (see gutMesh), the third input 
%             corresponds to which component of the velocity (longitudinal 
%             or transverse) and the last input corresponds to which frame.
%         - thetas: A 2xM array of slopes corresponding to the M columns.
%             Each slope is calculated as the perpedicular to the 
%             centerline.

function [gutMesh, mSlopes, gutMeshVels, gutMeshVelsPCoords, thetas] = interpolatePIVVectorsInMask(curDir, expDir, imageType, resReduce, rawPIVOutputName, maskFileOutputName)

%% Initialize variables
splineNFineness=10000;

%% Initialize mesh variables
% Load data
load(strcat(curDir, filesep, rawPIVOutputName,'_Current.mat'));
load(strcat(curDir, filesep, maskFileOutputName, '_Current.mat'));
ex=x{1}; %#ok as it should be loaded, represents the x positions of the gutMesh
why=y{1}; %#ok as it should be loaded, represents the x positions of the gutMesh

% Define NU and NV, the number of points anterior/posterior along the gut
% and the number of points dorsal/ventral to the gut
psIn=inpolygon(ex(:),why(:),gutOutlinePoly(:,1),gutOutlinePoly(:,2)); %#ok as it should be loaded
logicPsIn = reshape(psIn,size(ex));
NVDist=sum(logicPsIn,1); % Assumes gut is roughly horizontal
NV=2*ceil(mean(NVDist)/2);
NUDist=sum(logicPsIn,1);
NUDist(NUDist>0)=1;
NU=sum(NUDist); % "Raw" NU, but I want it nicer
NUrem=idivide(uint8(NU),10);
NU=double((NUrem+1)*10); % Just because I want an even, divisible by 10, NV

% Initialize the gutMesh variables (the variables containing spatial info)
finalExesTop=zeros(1,NU+1);
finalExesBottom=zeros(1,NU+1);
gutMesh=zeros(NV,NU,2);
mSlopes=zeros(NU,2,2); % of the form (position down gut, dx or dy, top or bottom)

% Interpolation won't work without removing duplicates (gutMiddlePolyTop
% and gutMiddlePolyBottom should be the same; historical reason for two
% variables comes from two middle regions being defined)
gutMiddlePolyUnTop=unique(gutMiddlePolyTop,'rows');
gutMiddlePolyUnBottom=unique(gutMiddlePolyBottom,'rows');

% Interpolate data, but only for each column
dUTop=ceil(size(gutMiddlePolyUnTop,1)/(NU+1));
dUBottom=ceil(size(gutMiddlePolyUnBottom,1)/(NU+1));
exTop=gutMiddlePolyUnTop(1:dUTop:end,1);
exBottom=gutMiddlePolyUnBottom(1:dUBottom:end,1);
whyTop=gutMiddlePolyUnTop(1:dUTop:end,2);
whyBottom=gutMiddlePolyUnBottom(1:dUBottom:end,2);
csTop=spline(exTop,whyTop);
csBottom=spline(exBottom,whyBottom);

% Subdivide the spline into NU even sections (so NU+1 vertices)
dExTop=(exTop(end)-exTop(1))/splineNFineness;
dExBottom=(exBottom(end)-exBottom(1))/splineNFineness;
exesTop=exTop(1):dExTop:exTop(end);
exesBottom=exBottom(1):dExBottom:exBottom(end);
curvezTop=ppval(csTop,exesTop);
curvezBottom=ppval(csBottom,exesBottom);
STotalTop=sum(sqrt(dExTop^2+diff(curvezTop).^2)); % Obtain the total length of the curve drawn
STotalBottom=sum(sqrt(dExBottom^2+diff(curvezBottom).^2)); % Obtain the total length of the curve drawn
dSTop=STotalTop/NU; % Arc length of N subdivisions
dSBottom=STotalBottom/NU; % Arc length of N subdivisions
partialSSumTop=cumsum(sqrt(dExTop^2+diff(curvezTop).^2)); % Partial sum of arcs
partialSSumBottom=cumsum(sqrt(dExBottom^2+diff(curvezBottom).^2)); % Partial sum of arcs
secSNumTop=uint32(floor(partialSSumTop/dSTop)); % Evenly number exes by arc length
secSNumBottom=uint32(floor(partialSSumBottom/dSBottom)); % Evenly number exes by arc length
for i=1:NU+1
    tempsNTop=find(secSNumTop==i-1);
    tempsNBottom=find(secSNumBottom==i-1);
    finalExesTop(i)=exesTop(tempsNTop(1)); % Vector of the y's in exes for which the y's in secSNum are first unique
    finalExesBottom(i)=exesBottom(tempsNBottom(1));
end

fullPolyX=[gutOutlinePoly(:,1); gutOutlinePoly(1,1)];
fullPolyY=[gutOutlinePoly(:,2); gutOutlinePoly(1,2)];
finalWhysTop=ppval(csTop,finalExesTop);
finalWhysBottom=ppval(csBottom,finalExesBottom);

% Load a representational image
baseFilenames = dir(strcat(expDir, filesep, imageType));
info = imfinfo([expDir filesep baseFilenames(1).name]);
y1=1;
y2 = info(1).Height/resReduce;

%% Generate mesh based on orthogonal vectors
for i=1:NU
    
    % Obtain mathematical function of the line from which centers will come
    curOrthMTop=-(finalExesTop(i+1)-finalExesTop(i))/(finalWhysTop(i+1)-finalWhysTop(i));
    curOrthMBottom=-(finalExesBottom(i+1)-finalExesBottom(i))/(finalWhysBottom(i+1)-finalWhysBottom(i));% Orthogonal slope
    midPointTop=[(finalExesTop(i+1)+finalExesTop(i))/2,(finalWhysTop(i+1)+finalWhysTop(i))/2 ];
    midPointBottom=[(finalExesBottom(i+1)+finalExesBottom(i))/2,(finalWhysBottom(i+1)+finalWhysBottom(i))/2 ];
    bTop=-curOrthMTop*midPointTop(1)+midPointTop(2);
    bBottom=-curOrthMBottom*midPointBottom(1)+midPointBottom(2);
    x1Top=(y1-bTop)/curOrthMTop;
    x1Bottom=(y1-bBottom)/curOrthMBottom;
    x2Top=(y2-bTop)/curOrthMTop;
    x2Bottom=(y2-bBottom)/curOrthMBottom;
    
    % Find their intersection with the gut edge previously drawn
    [exIntTop,whyIntTop]=polyxpoly([x1Top, x2Top],[y1, y2],fullPolyX,fullPolyY);
    [exIntBottom,whyIntBottom]=polyxpoly([x1Bottom, x2Bottom],[y1, y2],fullPolyX,fullPolyY);
    % This looks complicated: Keep in mind the intersection gives two
    % points, so bottom/top may refer to which gut midline or it may refer
    % to which intersection is in question. Further obfuscated by images
    % being upside down...
    topYIBottom=max(whyIntBottom);
    topCorXIBottom=exIntBottom(whyIntBottom==topYIBottom);
    bottomYITop=min(whyIntTop);
    bottomCorXITop=exIntTop(whyIntTop==bottomYITop);
    dYIBottom=topYIBottom-midPointBottom(2);
    dXIBottom=topCorXIBottom-midPointBottom(1);
    dYITop=bottomYITop-midPointTop(2);
    dXITop=bottomCorXITop-midPointTop(1);
    mSlopes(i,1,1)=dXITop;
    mSlopes(i,1,2)=-dXIBottom; % minus for inconsistency
    mSlopes(i,2,1)=dYITop;
    mSlopes(i,2,2)=-dYIBottom; % minus for inconsistency
    
    % Designate mesh locations (divide in to upper and lower for ease of indexing
    for j=1:NV/2
        
        % Top mesh locations
        gutMesh(NV-j+1,i,1)=round(midPointTop(1)+j*dXITop/(NV/2+1)); % curOrthM ratio for sign
        gutMesh(NV-j+1,i,2)=round(midPointTop(2)+j*dYITop/(NV/2+1));
                
        % Bottom mesh locations
        gutMesh(j,i,1)=round(midPointBottom(1)+j*dXIBottom/(NV/2+1)); % curOrthM ratio for sign
        gutMesh(j,i,2)=round(midPointBottom(2)+j*dYIBottom/(NV/2+1));
        
    end
    
end

% In the event there was no postprocessing?
if(~exist('u_filt','var'))
    u_filt=u;
end
if(~exist('v_filt','var'))
    v_filt=v;
end


%% Interpolate velocities
nT=size(u_filt,1)-1; % For nT frames, there should have only been nT-1 elements, instead they made the nTth element [] for some reason...
gutMeshVels=zeros(size(gutMesh,1),size(gutMesh,2),size(gutMesh,3),nT); % Ordered U=1, V=2

% Progress bar
progtitle = sprintf('Interpolatio');
progbar = waitbar(0, progtitle);  % will display progress
prevPercent = uint8(0);

for i=1:nT 
    
    if(uint8(i/nT*100)>prevPercent)
        prevPercent = uint8(i/nT*100);
        % Progress bar update
        waitbar(double(prevPercent)/100, progbar, ...
            strcat(progtitle, sprintf('n %d%% done', prevPercent)));
    end
    
    Xq=gutMesh(:,:,1);
    Yq=gutMesh(:,:,2);
    ex=x{i};
    why=y{i};
    V=v_filt{i};
    U=u_filt{i};
    SIU=scatteredInterpolant(ex(:),why(:),U(:));
    SIV=scatteredInterpolant(ex(:),why(:),V(:));
    mU=SIU(Xq(:),Yq(:));
    mV=SIV(Xq(:),Yq(:));
    gutMeshVels(:,:,1,i)=reshape(mU,size(Xq));
    gutMeshVels(:,:,2,i)=reshape(mV,size(Xq));
    
end

close(progbar);

%% Map to local coordinates
gutMeshVelsPCoords=gutMeshVels;
thetaStarTop=atan2(squeeze(mSlopes(:,2,1)),squeeze(mSlopes(:,1,1)));
thetaStarBottom=atan2(squeeze(mSlopes(:,2,2)),squeeze(mSlopes(:,1,2)));
thetaTop=thetaStarTop-pi/2;
thetaBottom=thetaStarBottom-pi/2;
thetas=[thetaTop, thetaBottom];

for i=1:size(thetaTop,1)
    gutMeshVelsPCoords(1:end/2,i,1,:)=gutMeshVels(1:end/2,i,1,:)*cos(thetaTop(i))+gutMeshVels(1:end/2,i,2,:)*sin(thetaTop(i));
    gutMeshVelsPCoords((end/2+1):end,i,1,:)=gutMeshVels((end/2+1):end,i,1,:)*cos(thetaBottom(i))+gutMeshVels((end/2+1):end,i,2,:)*sin(thetaBottom(i));
    gutMeshVelsPCoords(1:end/2,i,2,:)=-gutMeshVels(1:end/2,i,1,:)*sin(thetaTop(i))+gutMeshVels(1:end/2,i,2,:)*cos(thetaTop(i));
    gutMeshVelsPCoords((end/2+1):end,i,2,:)=-gutMeshVels((end/2+1):end,i,1,:)*sin(thetaBottom(i))+gutMeshVels((end/2+1):end,i,2,:)*cos(thetaBottom(i));
end

end