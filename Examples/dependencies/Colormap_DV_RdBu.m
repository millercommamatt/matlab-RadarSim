function [cm] = Colormap_DV_RdBu(nc, zeroPoint)
%Colormap_DV_RdBu This function generates a red-to-light gray-to-blue
%colormap that is suitable for doppler velocity plots or any situation
%where a diverging colormap is required. The colormap is based on the
%Brewer RdBu diverging colormap.
%
%INPUTS (all optional)
%  nc = number of color levels; default = 64
%  zeroPoint = The location of the zero or the light gray, central color.
%  If this value is between 0 and 1, it's treated as the fractional
%  location along the color map. If the value is an integer greater than
%  one and less than the total number of colors, it's interpreted as the
%  index where the light gray value should be set.
%
%OUTPUTS
%  cm = a colormap [nc x 3 double array]
%
%References:
%   http://mkweb.bcgsc.ca/brewer/swatches/brewer-palettes-swatches.pdf
%   http://colorbrewer2.org
%
%Written by: Matthew Miller, 2018

%Input validation
if nargin<1
    nc=64;
elseif nc < 3
    error('too few colors specified\n')
end

RdBu = ...
    [0.403921568627451,0,0.121568627450980;...
    0.698039215686275,0.0941176470588235,0.168627450980392;...
    0.792156862745098,0,0.125490196078431;
    0.839215686274510,0.376470588235294,0.301960784313725;...
    0.937254901960784,0.541176470588235,0.384313725490196;
    0.956862745098039,0.647058823529412,0.509803921568627;...
    0.992156862745098,0.858823529411765,0.780392156862745;...
    0.968627450980392,0.968627450980392,0.968627450980392;...
    0.819607843137255,0.898039215686275,0.941176470588235;...
    0.572549019607843,0.772549019607843,0.870588235294118;...
    0.403921568627451,0.662745098039216,0.811764705882353;
    0.262745098039216,0.576470588235294,0.764705882352941;...
    0.0196078431372549,0.443137254901961,0.690196078431373;
    0.129411764705882,0.400000000000000,0.674509803921569;...
    0.0196078431372549,0.188235294117647,0.380392156862745];

if exist('zeroPoint','var')
    if zeroPoint >0 && zeroPoint < 1 && nc >= 15
        zero_val = round(zeroPoint.*nc);
        cm1 = interp1(1:8,RdBu(1:8,:),linspace(1,8,zero_val),'pchip');
        cm2 = interp1(8:15,RdBu(8:15,:),linspace(8,15,nc-zero_val+1),'pchip');
    elseif isinteger(zeroPoint) && zeroPoint < nc
        cm1 = interp1(1:8,RdBu(1:8,:),linspace(1,8,zeroPoint),'pchip');
        cm2 = interp1(8:15,RdBu(8:15,:),linspace(8,15,nc-zeroPoint+1),'pchip');
    else
        error('Had trouble parsing the zeroPoint input.\n')
    end
    cm = [cm1;cm2(2:end,:)];
else
    switch nc
        case 3
            cm = RdBu([5 8 11],:);
        case 4
            cm = RdBu([3 6 10 13],:);
        case 5
            cm = RdBu([3 6 8 10 13],:);
        case 6
            cm = RdBu([2 5 7 9 11 14],:);
        case 7
            cm = RdBu([2 5 7 8 9 11 14],:);
        case 8
            cm = RdBu([2 4 6 7 9 10 12 14],:);
        case 9
            cm = RdBu([2 4 6 7 8 9 10 12 14],:);
        case 10
            cm = RdBu([1 2 4 6 7 9 10 12 14 15],:);
        case 11
            cm = RdBu([1 2 4 6 7 8 9 10 12 14 15],:);
        case 12
            cm = RdBu([1 2 3 4 6 7 9 10 12 13 14 15],:);
        case 13
            cm = RdBu([1 2 3 4 6 7 8 9 10 12 13 14 15],:);
        case 14
            cm = RdBu([1:7 9:15],:);
        case 15
            cm = RdBu;
        otherwise
            cm = interp1(1:15,RdBu,linspace(1,15,nc),'pchip');
    end
end


end

% %%%%%
% figure;plot_colorbar(10,'h','1 2 3',Colormap_DV_RdBu(64)); % red - blue
% figure;plot_colorbar(10,'h','3 2 1',fliplr(Colormap_DV_RdBu(64))); % blue - brown
% figure;plot_colorbar(10,'h','1 3 2',circshift(fliplr(Colormap_DV_RdBu(64)),1,2)); % red/pink - teal/green
% figure;plot_colorbar(10,'h','2 1 3',circshift(fliplr(Colormap_DV_RdBu(64)),2,2)); % green - purple
% figure;plot_colorbar(10,'h','3 1 2',circshift(Colormap_DV_RdBu(64),1,2));% green - pink/magenta
% figure;plot_colorbar(10,'h','2 3 1',circshift(Colormap_DV_RdBu(64),2,2));% blue/purple - mossy green
