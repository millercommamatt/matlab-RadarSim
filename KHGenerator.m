function [U, V, W] = KHGenerator(modelX,modelY,modelZ,waveLength,waveAltitude,waveThickness,waveDirection,velMagnitude,BGWindRange,phaseOffset)
%KHGenerator Generates 3D wind components for a fully planar
%Kelvin-Helmholtz wave.
%   KHGenerator(modelX,modelY,modelZ,waveLength,waveAltitude,waveThickness,waveDirection,velMagnitude,phaseOffset)
%   creates arrays of the 3D wind components of a planar Kelvin-Helmholtz
%   wave with user-specified characteristics for a specified 3D rectilinear
%   grid. The wave thickness is set using a normal distribution werein the
%   thickness is +/- 2 sigma. The background wind field concentrates shear
%   across the wave zone. The user defines the above and below wave values
%   in m/s and the vertical wind field varies according to a logistic
%   function. Adjusting the phase offset on successive calls can be used to
%   simulate wave propagation.
%
%   Input Arguments:
%   modelX          meshgrid array of X-coordinates in degrees longitude
%                   of the 3D wind values are calculated for
%
%   modelY          meshgrid array of X-coordinates in degrees latitude
%                   the 3D wind values are calculated for
%
%   modelX          meshgrid array of X-coordinates (down-up) in meters
%                   the 3D wind values are calculated for
%
%   waveLength      wavelength of the wave in meters
%
%   waveAltitude    height of the wave center in meters
%
%   waveThickness   vertical thickness of the wave in meters
%
%   waveDirection   direction of propagation of the wave in degrees
%                   clockwise from north
%
%   velMagnitude    maximum horizontal velocity value of the wave in m/s
%
%   BGWindRange     {optional) two element array with the below wave and
%                   above wave wind values in m/s.
%                   e.g. [below_wave above_wave] use [0 0] for no
%                   background wind field.
%
%   phaseOffset     (optional) phase offset of the wave in meters
%
%   Output Variables: 
%   U/V/W           meshgrid like arrays of the U, V, and W components of
%                   the 3D wind
%
%   Example:
%   Generate the wind components of a KH for the domain from -100 - -98
%   degrees longitude, 30 - 32 degrees latitude, and from 0 - 10,000 meters
%   altitude with 100 steps in the horizontal dimensions and 30 in the
%   vertical dimension. The wavelength will be 3,000 meters, the altitude
%   will be 2,500 meters, thickness will be 800 meters, the maximum
%   horizontal windspeed of the wave will be 4.2 m/s, the wave will be
%   oriented towards the west (270 deg.), the horizontal background wind
%   will be -2 m/s below the wave and 3 m/s above the wave, and wave phase
%   will not be offset (0 meters).
%
%   %Create spatial domain
%   [modelX, modelY, modelZ] = meshgrid(...
%      linspace(-100, -98, 100),...
%      linspace(30, 32, 100),...
%      linspace(0, 10000, 30));
%   %generate waves
%   [U,V,W] = KHGenerator(...
%      modelX, modelY, modelZ,...
%      3000, 2500, 800, 270, 4.2, [-2 3], 0);
%      
%
%   Written By: Dr. Matthew A. Miller, PhD, 2021
%
%   References: Grasmick, Coltin, and Bart Geerts. " Detailed Dual-Doppler 
%                 Structure of Kelvin–Helmholtz Waves from an Airborne 
%                 Profiling Radar over Complex Terrain. Part I: Dynamic 
%                 Structure", Journal of the Atmospheric Sciences 77, 5 
%                 (2020): 1761-1782, accessed May 7, 2021, 
%                 https://doi.org/10.1175/JAS-D-19-0108.1


% INPUT VALIDATION
if nargin < 8
    BGWindRange = [-velMagnitude velMagnitude];
end
if nargin < 10
    phaseOffset = 0;
end
validateattributes(waveLength, {'numeric'},{'scalar','>',0},mfilename,'Wavelength in meters',4)
validateattributes(waveAltitude, {'numeric'},{'scalar','>',0,'<=',50000},mfilename,'Wave Altitude',5)
validateattributes(waveThickness, {'numeric'},{'scalar','>',0,'<=',20000},mfilename,'Wave Thickness',6)
validateattributes(waveDirection, {'numeric'},{'scalar','>=',0,'<=',360},mfilename,'Wave Direction',7)

[ySize,xSize,zSize]=size(modelZ);

%compute X direction distance of the model domain at the box middle
%There is an expectation here for data that's lat/lon aligned. Model data
%may not meet that criteria.
xModelDist = distance(modelY(round(end/2),1,1),modelX(1,1,1),...
                      modelY(round(end/2),1,1), modelX(1,end,1),...
                      referenceEllipsoid('earth', 'm'));
yModelDist = distance(modelY(1,1,1),modelX(1,round(end/2),1),...
                      modelY(end,1,1), modelX(1,round(end/2),1),...
                      referenceEllipsoid('earth', 'm'));

numWaves = xModelDist/waveLength;

[X,Y,~] = meshgrid(linspace(0,2*pi,xSize),...
                   linspace(0,2*pi*(yModelDist/xModelDist),ySize),...
                   1:zSize);

Dir = compass2math(waveDirection);

% Make the wave

% Make the height mask which is a normal distribution shape with a max
% value of 1. The standard deviation is the waveThickness_parameter and is
% used to set the thickess of the wave as +/- 2 SD. 
waveThickness_parameter = waveThickness / (2 * 2); % 2 sigma
heightMask = normpdf(modelZ,waveAltitude,waveThickness_parameter)./(1/(waveThickness_parameter*sqrt(2*pi)));

% The horizontal and vertical waves are defined as sin/cos waves a quarter
% wavelength offset that then shift a quarter wavelength over the thickess
% of the wave from top to bottom per the schematic from Figure 11 of
% Grasmick and Geerts, 2020, JAS. The hight mask is used to set the
% perturbation towards zero since the functions used to define the
% perturbation would extend the full Z depth.
vert_phase_shift = (pi/2)/waveThickness;  % 1/4 wavelength over wave thickness
H_wave = heightMask.*velMagnitude.*sin(numWaves.*(cosd(Dir).*X + sind(Dir).*Y) + vert_phase_shift.*modelZ -(phaseOffset/waveLength)*2*pi);
W = heightMask.*velMagnitude.*cos(numWaves.*(cosd(Dir).*X + sind(Dir).*Y) - vert_phase_shift.*modelZ -(phaseOffset/waveLength)*2*pi);

% Define a background horizontal wind. The background wind is defined with
% components parallel to the shear vector (i.e. waveDirection). The user
% define the magnitudes above and below the wave. The magnitudes shift from
% the lower value to the upper value according to a logistic function
% across the depth of the wave.
BG_range = BGWindRange(2)-BGWindRange(1);
LogisticSteepness = 10 / waveThickness;
H_bg = (BG_range./(1+exp(-LogisticSteepness.*(modelZ-waveAltitude)))) + BGWindRange(1);

U = (H_bg + H_wave).*cosd(Dir);
V = (H_bg + H_wave).*sind(Dir);

end

