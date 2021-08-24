function [U,V,W] = waveGenerator(modelX,modelY,modelZ,waveLength,waveHeight,waveDirection,velMagnitude,phaseOffset,waveTilt)
%WAVEGENERATOR Generates the 3D wind components for a monochromatic plane
%wave.
%   WAVEGENERATOR(modelX,modelY,modelZ,waveLength,waveHeight,waveDirection,velMagnitude,phaseOffset,waveTilt)
%   creates arrays of the 3D wind components of a monochromatic plane wave
%   with user-specified characteristics for a specified 3D rectilinear
%   grid. The wave extends from the surface to the specified wave height.
%   The wave decays with height starting from the user-given height through
%   a depth of 10% of the total wave depth according to a reverse logistic
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
%   waveHeight      height of the wave top in meters
%
%   waveDirection   direction of propagation of the wave in degrees
%                   clockwise from north
%
%   velMagnitude    maximum horizontal velocity value of the wave in m/s
%
%   phaseOffset     (optional) phase offset of the wave in meters
%
%   waveTilt        (optional, phaseOffset required) vertical tilt of the
%                   wave in degrees from vertical
%
%   Output Variables: 
%   U/V/W           meshgrid like arrays of the U, V, and W components of
%                   the 3D wind
%
%   Example:
%   Generate the wind components of a planer wave for the domain from -100
%   - -98 degrees longitude, 30 - 32 degrees latitude, and from 0 - 10,000
%   meters altitude with 100 steps in the horizontal dimensions and 30 in
%   the vertical dimension. The wavelength will be 30,000 meters, the
%   height will be 2,500 meters, the maximum horizontal windspeed of the
%   wave will be 4.2 m/s, the wave will be oriented towards the west (270
%   deg.), and wave phase will not be offset (0 meters), and the wave tilt
%   will be 45 degrees (1:1 slope).
%
%   %Create spatial domain
%   [modelX, modelY, modelZ] = meshgrid(...
%      linspace(-100, -98, 100),...
%      linspace(30, 32, 100),...
%      linspace(0, 10000, 30));
%   %generate waves
%   [U,V,W] = waveGenerator(...
%      modelX, modelY, modelZ,...
%      30000, 2500, 270, 4.2, 0, 45);
%      
%
%   Written By: Dr. Matthew A. Miller, PhD, 2017

% INPUT VALIDATION
if nargin < 8
    phaseOffset = 0;
end
if nargin < 9
    waveTilt = 0;
end
validateattributes(waveLength, {'numeric'},{'scalar','>',0},mfilename,'Wavelength in meters',4)
validateattributes(waveHeight, {'numeric'},{'scalar','>',0,'<=',100000},mfilename,'Wave Height',5)
validateattributes(waveDirection, {'numeric'},{'scalar','>=',0,'<=',360},mfilename,'Wave Direction',6)
validateattributes(waveTilt, {'numeric'},{'scalar','>',-90,'<',90},mfilename,'Wave Tilt',9)


[ySize,xSize,zSize]=size(modelZ);

%compute X direction distance of the model domain at the box middle
%There is an expectation here for data that's lat/lon aligned. Model data
%may not meet that criteria.
xModelDist = distance(modelY(round(end/2),1,1),modelX(1,1,1), modelY(round(end/2),1,1), modelX(1,end,1), referenceEllipsoid('earth', 'm'));
yModelDist = distance(modelY(1,1,1),modelX(1,round(end/2),1), modelY(end,1,1), modelX(1,round(end/2),1), referenceEllipsoid('earth', 'm'));

xModelSpacing = xModelDist./xSize;
yModelSpacing = yModelDist./ySize;
zModelSpacing = modelZ(round(end/2),round(end/2),2) - modelZ(round(end/2),round(end/2),1);
%this will work badly when trying to fit model data because it's neither
%spaced evenly with height nor consistant from column to column

numWaves = xModelDist/waveLength;

[X,Y,Z] = meshgrid(linspace(0,2*pi,xSize),linspace(0,2*pi*(yModelDist/xModelDist),ySize),linspace(0,2*pi*(range(modelZ(:))/xModelDist)*numWaves,zSize));

Dir = compass2math(waveDirection);
hVel = velMagnitude.*sin(numWaves.*(cosd(Dir).*X + sind(Dir).*Y) - tand(waveTilt).*Z -(phaseOffset/waveLength)*2*pi);

Ubase = hVel.*cosd(Dir);
Vbase = hVel.*sind(Dir);

%decay wave with height using logistic function
decayDepth = waveHeight/10;%200;%in meters
decayVector = 1-(1./(1+exp(-(4/decayDepth).*(modelZ(1,1,:)-waveHeight))));
decayVector = reshape(decayVector,1,1,[]);
decayArray = repmat(decayVector,[ySize,xSize,1]);

U = decayArray.*Ubase;
V = decayArray.*Vbase;

%use pseudo mass continuity to create W component;
W = zeros(ySize,xSize,zSize);
%%% VERTICAL MOTION OMITTED UNTIL A GOOD SOLUTION CAN BE FOUND
%
% % Approximation of standard atmosphere density as a function of height in
% % meters. This is used to the pesudo mass continuity calculation.
% atmoDens = @(x) 3.307899E-9.*x.^2 - 1.136965E-4.*x + 1.221873;
% 
% for ii = 1:zSize
%     
%     if ii == 1
%         %W(:,:,ii) = (uFlux(:,:,ii).*(yModelSpacing*zModelSpacing) + vFlux(:,:,ii).*(xModelSpacing*zModelSpacing))./(xModelSpacing*yModelSpacing);
%         W(:,:,ii) = (uFlux(:,:,ii).*(yModelSpacing*zModelSpacing) + vFlux(:,:,ii).*(xModelSpacing*zModelSpacing)).*atmoDens(modelZ(:,:,ii))./(xModelSpacing*yModelSpacing.*atmoDens(modelZ(:,:,ii)+zModelSpacing.*0.5));
%     else
%         %W(:,:,ii) = ((uFlux(:,:,ii).*(yModelSpacing*zModelSpacing) + vFlux(:,:,ii).*(xModelSpacing*zModelSpacing) + W(:,:,ii-1).*(xModelSpacing*yModelSpacing))./(xModelSpacing*yModelSpacing)).*(1-(0.0091*zModelSpacing));%rough approximation of change in mass with height based on a standard atmosphere
%         W(:,:,ii) = (((uFlux(:,:,ii).*(yModelSpacing*zModelSpacing) + vFlux(:,:,ii).*(xModelSpacing*zModelSpacing)).*atmoDens(modelZ(:,:,ii))) + W(:,:,ii-1).*(xModelSpacing*yModelSpacing).*atmoDens(modelZ(:,:,ii)-zModelSpacing.*0.5))./(xModelSpacing*yModelSpacing.*atmoDens(modelZ(:,:,ii)+zModelSpacing.*0.5));
%         %W(:,:,ii) = ((uFlux(:,:,ii).*(yModelSpacing*zModelSpacing) + vFlux(:,:,ii).*(xModelSpacing*zModelSpacing) + W(:,:,ii-1).*(xModelSpacing*yModelSpacing))./(xModelSpacing*yModelSpacing)).*( atmoDens( modelZ(:,:,ii) ) ./ atmoDens( modelZ(:,:,ii-1) ) );
%     end
% end


