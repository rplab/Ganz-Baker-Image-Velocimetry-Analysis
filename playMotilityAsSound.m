% Converts motility into sound
% Assumes gutMesh, gutMeshVelsPCoords exist
%
% Variables: totalTimeFraction, Use 1 if all
%            fractionOfTimeStart, Use size(gutMeshVelsPCoords,4) if all
%            markerNumStart, Use 1 for anterior position
%            markerNumEnd, Use size(gutMesh,2) if all

function theSound = playMotilityAsSound(gutMeshVelsPCoords, totalTimeFraction, fractionOfTimeStart, markerNumStart, markerNumEnd)

% Initialize variables
sinDensity = 100;
ordinateValues= int16((size(gutMeshVelsPCoords,4)/fractionOfTimeStart):(size(gutMeshVelsPCoords,4)/(fractionOfTimeStart)+size(gutMeshVelsPCoords,4)/totalTimeFraction-1));
surfaceValues=squeeze(-mean(gutMeshVelsPCoords(:,markerNumStart:markerNumEnd,1,ordinateValues),1));
theSound = zeros(1,sinDensity*size(surfaceValues, 2));

% Generate sound
for i=1:size(surfaceValues, 2)
    for j=1:size(surfaceValues, 1)
        sinSound = sin(linspace(0, 3.141592*j/2, sinDensity));
        theSound((sinDensity*(i-1) + 1):sinDensity*i) = theSound((sinDensity*(i-1) + 1):sinDensity*i) + sinSound*surfaceValues(j,i);
    end
end

end