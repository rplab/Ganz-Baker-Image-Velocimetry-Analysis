# Summary of Toolbox and function dependencies

Summary of Toolbox and function dependencies of analyzeMotility.m and its sub-functions, assessed using the excellent "fdep_21jun2010" utility (File Exchange. Author?) See also my "script_to_list_toolboxes.m")

Raghuveer Parthasarathy
August 11, 2018

### Toolboxes Used:

- Mapping Toolbox
- Signal Processing Toolbox
- Image Processing Toolbox

### Image Processing Toolbox

Lots of functions use parts of this; I'm not going to list them, since I expect that this toolkit is necessary.

### Signal Processing Toolbox

- gutFreqWaveSpeedFinder.m -- uses designfilt.p, filtfilt.m
- obtainMotilityParameters.m -- uses designfilt.p, filtfilt.m, xcorr.m

### Mapping Toolbox [obsolete]

- interpolatePIVVectorsInMask.m . Uses polyxpoly.m

Note that polyxpoly.m calls checkxy.m, and intersectLineSegments.m, so these are also needed.

See, however, the next item.

### Avoiding the Mapping toolbox

To avoid needing polyxpoly.m, etc., get "intersections.m" from the FileExchange:
https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections  (Douglas Schwarz)

Called in interpolatePIVVectorsInMask.m

**Path**

e.g. (on my computer) 
'/Applications/MATLAB_R2015a.app/toolbox/signal/signal/designfilt.p'.