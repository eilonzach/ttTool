function [ZZ] = gauss1D(X,sX,xx)
% function [ZZ] = gauss1D(X,sX,xx)
% 
% function to create a time series with a peak at X with x-width
% described by sX and resolved onto a line with dimension xx.


ZZ = exp(-( ((xx - X).^2)/(2*sX^2)));