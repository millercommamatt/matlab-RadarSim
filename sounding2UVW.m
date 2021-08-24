function [U, V, W] = sounding2UVW(u, v, numx, numy)
%SOUNDING2UVW Creates meshgrid arrays of U/V/W winds from an input
%sounding.
%
%   [U, V, W] = sounding2UVW(u, v, numx, numy) Constructs meshgrid format
%   arrays of U, V, and W winds from an input sounding u and v wind values
%   and a user-specified x and y grid spacing. The spacing in the W
%   direction is determined by the length of the u and v inputs.
%
%   Assume W = 0 for all points since sounding typically don't include
%   vertical velocity.
%
%   Input Arguments:
%   U         Vector of the U wind component. Length is equal to the number
%             of heights in the sounding.
%
%   V         Vector of the U wind component. Length is equal to the number
%             of heights in the sounding. The dimensions of U and V must
%             match.
%
%   NUMX      Positive, interger scalar of the number of desired 
%             coordinates in the x-direction.
%
%   NUMY      Positive, interger scalar of the number of desired 
%             coordinates in the y-direction.
%
%   Output Arguments:
%   U, V, W   numy x numx x length(U) array of U, V, and W wind components
%             respectively.
%
%   Examples: Load sounding data, convert that speed direction data to U
%             and V components and then create meshgrid-format U V W
%             arrays for 75 x-direction points and 150 y-direction points.
%
%      soundingData = load('sounding.m')
%      
%      soundingData =
%         windSpeed = [1 x 100 double]
%         windDir = [1 x 100 double]
%         height = [1 x 100 double]
%
%      soundingU = soundingData.windSpeed .* ...
%         cosd(compass2math(soundingData.windDir));
%      soundingV = soundingData.windSpeed .* ...
%         sind(compass2math(soundingData.windDir));
%
%      [U, V, W] = sounding2UVW(soundingU, soundingU, 75, 150);
%
%      size(radarU)
%
%      ans =
%         150         75        100
%
%   Written By: Dr. Matthew A. Miller, PhD, 24 May 2017

%VALIDATE INPUT ATTRIBUTES
validateattributes(u, {'numeric'},{'vector'},mfilename,'u wind',1)
validateattributes(v, {'numeric'},{'vector','size',size(u)},mfilename,'v wind',2)
validateattributes(numx, {'numeric'},{'scalar','integer','positive'},mfilename,'numx',3)
validateattributes(numy, {'numeric'},{'scalar','integer','positive'},mfilename,'numy',4)

% set permute order
if isrow(u)
    permuteOrder = [1 3 2];
else
    permuteOrder = [2 3 1];
end

U = repmat(permute(u,permuteOrder),[numy numx]);
V = repmat(permute(v,permuteOrder),[numy numx]);
W = zeros(numy,numx,length(u));

%check output
%size of a meshgrid array should be [numy, numx]

end

