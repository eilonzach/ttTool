function [ x,y ] = Lowhemis_conv( phi,theta,rmax )
% [ x,y ] = Lowhemis_conv( phi,theta,rmax )
%
% function to convert points described by angles phi (azimuth from N) and
% theta (angle from horizontal) on a lower hemisphere projection to x,y
% cartesian coordinates. E.g. a horizontal line going north would be at 
%     phi,theta = [90,90] and x,y = [0,1]   (or [0,rmax])
% The projection is an equal area projection (so 45 degree dip plots
% halfway between the origin and rmax
%
% INPUTS:
%   phi - azimuth from North/up the page (degrees)
%   theta - azimuth from horizontal, like dip (degrees)
%   rmax - maximum radius from the origin - for theta=0 (default is 1)
% 
% OUTPUTS:
%   x - horizontal cartesian coordinate (ordinate)
%   y - vertival cartesian coordinate (abcissa)

if nargin < 3
    rmax = 1;
end

r_ = rmax * sind((90-theta)./2)./sind(45);
x = sind(phi).*r_;
y = cosd(phi).*r_;

end

