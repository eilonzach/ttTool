function [ da ] = ang_diff( a1,a2,unit,mode )
%function [ da ] = ang_diff( a1,a2,unit,mode )
%   Function to work out the absolute angular difference between two angles
%   = a solution to the wraparound problem
%
% INPUTS
% a1 = the first angle (can be a vector of angles)
% a2 = the second angle
% unit = 'degrees' or 'radians' (default is degrees)
% mode = 'direct' or 'axial' (default is direct) treat data as linear or
% axial
% 
% OUTPUT
% da = the absolute angular difference between a1 and a2, in input units

if nargin<3
    unit='degrees';
end
if nargin<4
    mode='direct';
end

if strcmp(unit,'radians')==1
    a1 = (180/pi).*a1;
    a2 = (180/pi).*a2;
elseif strcmp(unit,'degrees')~=1
    error('Are these angles in degrees or radians?')
end

if strcmp(mode,'axial')
    a1 = mod(2.*a1,360)./2;
	a2 = mod(2.*a2,360)./2;
end

C = cosd(a1)+cosd(a2);
S = sind(a1)+sind(a2);
R=(C.^2 + S.^2);
da=180-acosd(1 - R/2);

if strcmp(mode,'axial')
    da1=da;
    a2=a2+180;
    C = cosd(a1)+cosd(a2);
    S = sind(a1)+sind(a2);
    R = (C.^2 + S.^2);
    da2=180-acosd(1 - R/2);
    if iscolumn(a1)
    da=min([da1,da2],[],2);
    else
    da=min([da1;da2],[],1);
    end
end

if strcmp(unit,'radians')
da = d2r(da);
end
    
end

