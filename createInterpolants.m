function [ interpolants ] = createInterpolants( modelX,modelY,modelZ,modelU,modelV,modelW)
%CALCULATE INTERPOLANTS creates the interpolation functions needed to
%calculate doppler velocity from a given model grid.
%   CREATEINTERPOLANTS(modelX,modelY,modelZ,modelU,modelV,modelW)
%   creates the interpolation functions needed to calculate doppler
%   velocity from a given model grid. The output structure is intended to
%   be used with the calculateDV function.
%
%   Input Arguments:
%   modelX/modelY/modelZ   3D meshgrid-like array of X/Y/Z coordinates of
%                          the field radar radial velocity will be 
%                          calculated for. X and Y dimensions are in 
%                          degrees longitufe and latitude respectively. The
%                          Z dimension inheight above geoid in meters.
%
%   modelU/modelV/modelW   3D meshgrid-like array of U/V/W windcomponents
%                          in meters per second.
%
%   Output Variables:
%   interpolants           Structure containing interpolation class objects
%                          for the U/V/W wind components. May be
%                          griddedinterpolant or scatteredinterpolant
%                          object depending on the nature of the input
%                          data.
%
%   See also meshgrid, ndgrid, griddedInterpolant, scatteredInterpolant
%
%   Written By: Dr. Matthew A. Miller, PhD, May 2017


%this is all for changing from meshgrid format to ndgrid format.
% The vertical coordinate is in meters while the x and y coordiantes are in
% degrees lat/lon. The typical range for lat/lon values are typically less
% than 10 while the typical range for Z can be 10000 or more. For the
% interpolation to work best, all the dimensions should have the same
% number of digits. So, the Z values are divided by 10000.
P = [2 1 3];
X = permute(modelX, P);
Y = permute(modelY, P);
Z = permute(modelZ, P)./10000;
U = permute(modelU, P);
V = permute(modelV, P);
W = permute(modelW, P);

%%%%%%%%%%%%%
%   TO DO   %
%%%%%%%%%%%%%
%
% Interpolation could happen in meter space. That means the lat and
% lon values need to be converted to meter offset a la aer2ecef.
%
%%%%%%%%%%%%%
% END TO DO %
%%%%%%%%%%%%%

%Inputs are assumed to be gridded. If that fails, they'll be treated as
%scattered
try
    interpolants.u = griddedInterpolant(X,Y,Z,U,'linear','none');
    interpolants.v = griddedInterpolant(X,Y,Z,V,'linear','none');
    interpolants.w = griddedInterpolant(X,Y,Z,W,'linear','none');
catch
    interpolants.u = scatteredInterpolant(X(:),Y(:),Z(:),U(:),'linear','none');
    interpolants.v = interpolants.u;
    interpolants.w = interpolants.u;
    interpolants.v.Values = V(:);
    interpolants.w.Values = W(:);
end

end

