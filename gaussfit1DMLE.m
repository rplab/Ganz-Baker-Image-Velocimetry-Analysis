% gaussfit1DMLE.m
%
% function to fit a 1D Gaussian using Maximum-Likelihood
% Estimation (MLE):  z(x) = A*exp(-((x-x0)^2)/2/sigma^2) + offset;
% Fits to a Gaussian with a constant (background) offset and
% *Poisson-distributed* noise.
%
% Checks if fit center is outside of the range; if so, output x0=centroid.
% (MLE can give nonsensical results, esp. at low signal/noise)
%
% Input: 
%    z : 1D array 
%    tolerance in z for fitting (default 1e-6 if empty)  
%    params0 : [optional] starting values of the parameters
%              If not input or empty, default values
%              1 - constant offset (minimal value of z)
%              2 - Gaussian center (x-coordinate of the "center of mass")
%              3 - sigma   (std. dev. of x range)
%              4 - amplitude (max. z - min. z)
%    px : x values (optional); default is 1:length(z). If these are input, params0
%        must be input, since default calculation won't work.
% Output:
%    A  : Gaussian amplitude
%    x0 : Gaussian center .
%    sigma: Std dev. of Gaussian (assumed same for x)
%    offset : constant offset
%
% Raghuveer Parthasarathy 
% August 3, 2017; based on gaussfit2DMLE.m
% October 24, 2021: major fix to optimization function, mistakenly included 
%    the sum that was used in 2D and not applicable in 1D
% Last modified July 12, 2023

function [A, x0, sigma, offset] = gaussfit1DMLE(z, tolz, params0, px)

z = z(:);

nx = length(z);
if ~exist('px', 'var') || isempty(px)
    px = (1:nx)';
end

% defaults for initial parameter values, and lower and upperbounds
if ~exist('tolz', 'var') || isempty(tolz)
    tolz = 1e-6;
end
if ~exist('params0', 'var') || isempty(params0)
    params0 = [min(z(:)), sum(px(:).*z(:))/sum(z),std(px), max(z(:))-min(z(:))];
end

% More fitting options
fminoptions.TolFun = tolz;  %  % MATLAB default is 1e-6
fminoptions.TolX = 1e-6';  % default is 1e-6
fminoptions.Display = 'off'; % 'off' or 'final'; 'iter' for display at each iteration
fminoptions.LargeScale = 'off';  % use the medium-scale algorithm,
   % since I'm not supplying the gradient
   
% Add an offset, to avoid the chance of small numbers causing problems for
% the estimator function (only makes a difference for small signal/noise,
% less than about 6.)
zoff = 0.1;  % Does choice of zoff matter?
z = z+zoff;
params0(1) = params0(1) + zoff;

params = fminunc(@(P) objfun(P,px,z),params0,fminoptions);
    
A = params(4);
x0 = params(2);
sigma = abs(params(3));  % unconstrained minimization; sigma may be negative
offset = params(1) - zoff;

% Check if particle center is outside of image; return simple centroid if so
if (x0<min(px)) || (x0>max(px)) 
    x0 = sum(z.*px, 'all')/sum(z(:));
    % disp(px)
    % disp(z)
    % disp(z.*px)
    % disp(sum(z.*px))
    % figure; plot(px, z, 'kx-')
    % disp(x0)
    % disp('here 1')
    % disp('gaussfit1DMLE out of bounds: returning (x0) = centroid!')
end

end

    function negL = objfun(params,px,z)
    %  Log-likelihood for Poisson-distributed intensity values, with a
    %  Gaussian spatial distribution
        gaussprob = params(4)*exp(-(px(:) - params(2)).^2/2/params(3)/params(3)) + params(1);
        Lk = (z(:)).*log(gaussprob) - gaussprob;  % log-likelihood; Poisson-distr.
        negL = -sum(Lk(:));
    end
    