function [ radarStruct ] = createRadarStructure(radarLon, radarLat, beamWidth, tilts, range, gateSpacing, antennaHeight)%,overSampleFactor)
%CREATERADARSTRUCTRUE Create the radar structure object used in RadarSim.
%   CREATERADARSTRUCTRUE(radarLon, radarLat, beamWidth, tilts, range, gateSpacing, antennaHeight)
%   creates a structure that contains the georegistered coordinates of the
%   centers of all the radar gates and the beam-parallel unit vectors
%   necessary to calculate the doppler velocity. All input variables are
%   required. All units are in meters or degrees (compass relative).
%
%   Input Arguments:
%   RADARLON        Scalar value, in degrees, of the longitude of the radar
%                   location.
%
%   RADARLAT        Scalar value, in degrees, of the latitude of the radar
%                   location.
%
%   BEAMWIDTH       Scalar value, in degrees, of the half-power beam width.
%                   Values between 0.05 and 20 degrees supported.
%
%   TILTS           Numeric vector of increasing values, in degrees, of the
%                   elevation of the radar tilts. Values between -10 and 90
%                   degrees supported.
%
%   RANGE           Scalar value, in meters, of the maximum radar range.
%                   Values up to 500 Km (500000 m) supported.
%
%   GATESPACING     Scalar value, in meters, of the spacing between range
%                   gates. Values greater than 0 and less than half the
%                   range are supported.
%
%   ANTENNAHEIGHT   Scalar value, in meters, of the radar antenna above the
%                   geoid (wgs84). Values between -500 and 8000 supported.
%
%   Output Arguments:
%   RADARSTRUCT    A structure array containing the following fields:
%      
%      latRadar/lonRadar/hRadar - meshgrid format 3D arrays of the
%      latitude and longitude (in degrees) and the height (in meters) of
%      each radar sample volume (wgs84 reference). The array size is #tilts
%      x #rays x #gates.
%      
%      uUnitRadar/vUnitRadar/wUnitRadar - The U/V/W components of the
%      beam-parallel unit vector for each sample volume. These are required
%      to calculate the doppler velocity.
%
%   Example: Define the radar structure for a radar at 85.2 degrees W
%            longitude, 25 degrees N latitude, with a 0.9 degree beam width,
%            two tilts and 0.5 and 2.5 degrees, a maximum range of 250 Km,
%            a gate spacing of 250 meters, and an antenna height 124 meters
%            above the geoid.
%   
%      radarStruct = createRadarStructure(-85.2, 25, 0.9, [0.5 2.5],...
%         250000, 250, 124);
%
%   See also validatestring, aer2geodetic, sph2cart
%
%   Written By: Dr. Matthew A. Miller, PhD, May 2017

% ADD INPUT VALIDATION 
validateattributes(radarLon, {'numeric'},{'scalar','>=',-180,'<=',180},mfilename,'Radar Lon.',1)
validateattributes(radarLat, {'numeric'},{'scalar','>=',-90,'<=',90},mfilename,'Radar Lat.',2)
validateattributes(beamWidth, {'numeric'},{'scalar','>=',0.05,'<=',20},mfilename,'Beam width',3)
validateattributes(tilts, {'numeric'},{'vector','increasing','>=',-10,'<',90},mfilename,'Tilt(s)',4)
validateattributes(range, {'numeric'},{'scalar','>',-0,'<=',500000},mfilename,'Radar Range',5)
validateattributes(gateSpacing, {'numeric'},{'scalar','>',0,'<=',range/2},mfilename,'Radar Gate Spacing',6)
validateattributes(antennaHeight, {'numeric'},{'scalar','>',-500,'<=',8000},mfilename,'Radar Antenna Height (above geoid)',7)
%validateattributes(overSampleFactor, {'numeric'},{'scalar','>',0,'<=',10,'integer'},mfilename,'overSampleFactor',8)
    

%generate a radar x/y/z

azimuths = 0:beamWidth:360;%beamWidth:beamWidth:360;

%[lat,lon,h] = aer2geodetic(az,elev,slantRange,lat0,lon0,h0,spheroid);
[az,elev,slantRange]=meshgrid(azimuths,tilts,gateSpacing:gateSpacing:range);
%[latRadar,lonRadar,hRadar] = aer2geodetic(az,elev,slantRange,radarLoc(2),radarLoc(1),radarLoc(3),referenceEllipsoid('earth','meter'));
% this is handy but doesn't accout for refraction
%
%Reinhart uses the 4/3 earth radius appromization to account for both
%curvature and refraction. If you use aer2geodetic with a 4/3 earth radius
%oblate spheroid, the results match the formula given by Rinhart
earthSphere = referenceEllipsoid('earth','meter');
% refractSphere = oblateSpheroid;
% refractSphere.SemimajorAxis=earthSphere.SemimajorAxis.*(4/3);
% refractSphere.SemiminorAxis=earthSphere.SemiminorAxis.*(4/3);

%[latRadar,lonRadar,~] = aer2geodetic(az,elev,slantRange,radarLoc(2),radarLoc(1),antennaHeight),earthSphere);
%[~,~,hRadar] = aer2geodetic(az,elev,slantRange,radarLoc(2),radarLoc(1),radarLoc(3),refractSphere);
%[radarStruct.latRadar,radarStruct.lonRadar,radarStruct.hRadar] = aer2geodetic(az,elev,slantRange,radarLat,radarLon,antennaHeight,refractSphere);

%%%
localEarthRadius = rcurve('transverse',referenceEllipsoid('earth','meters'),radarLat,'degrees');
R = (4/3)*localEarthRadius;
radarStruct.hRadar = sqrt(slantRange.^2+R^2+2.*slantRange.*R.*sind(elev))-R+antennaHeight;%beam height above geoid in m
s = R .* asin((slantRange.*cosd(elev))./(R + radarStruct.hRadar));%Distance along geoid from radar in m
radarStruct.Xdist = s .* sind(az);%in meters
radarStruct.Ydist = s .* cosd(az);%in meters
[radarStruct.latRadar,radarStruct.lonRadar] = reckon(radarLat,radarLon,s./1000,az,referenceEllipsoid('earth','km'));
%%%


%Due to refraction, we need to compute the beam-relative azimuth and
%elevation in order to calculate the doppler velocity as the radar sees it.

[beamAz,beamElev,~] = geodetic2aer(radarStruct.latRadar(:,:,2:end),radarStruct.lonRadar(:,:,2:end),radarStruct.hRadar(:,:,2:end),...
    radarStruct.latRadar(:,:,1:end-1),radarStruct.lonRadar(:,:,1:end-1),radarStruct.hRadar(:,:,1:end-1),earthSphere);
beamElev = cat(3,elev(:,:,1),beamElev);
beamAz = cat(3,az(:,:,1),beamAz);

%compute radar unit vector
%[xUnitRadar,yUnitRadar,zUnitRadar] = sph2cart(deg2rad(compass2math(az)),deg2rad(beamElev),1);
[radarStruct.uUnitRadar,radarStruct.vUnitRadar,radarStruct.wUnitRadar] = sph2cart(deg2rad(compass2math(beamAz)),deg2rad(beamElev),1);


end

%LOCAL FUNCTION
function [ mathDeg ] = compass2math( compassDeg )
%COMPASS2MATH Converts compass degrees (0 at North) to math degrees(0 at
%"east")
%   Detailed explanation goes here

mathDeg = compassDeg;
mathDeg(compassDeg>=0 & compassDeg<90) = 90-compassDeg(compassDeg>=0 & compassDeg<90);
mathDeg(compassDeg>=90 & compassDeg<=360) = 450-compassDeg(compassDeg>=90 & compassDeg<=360);

end

