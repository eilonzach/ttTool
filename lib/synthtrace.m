function [ xx ] = synthtrace(samplength,impulselength,amp,dt,waveform,startfraction,rolloff)
%function [ xx ] = synthtrace(samplength,impulselength,amp,dt,startfraction,waveform,rolloff)
% 
% Function to output a synthetic waveform pulse; a time series with an
% artificial pulse in it
%
%INPUTS:
% samplength    - length (seconds) of whole time series
% impulselength - length (seconds) of pulse =(?) characteristic period
% amp           - amplitude of pulse - NB must adjust to get the right power 
% dt            - time between samples = 1/frequency, i.e. 1/samprate 
% startfraction - time fraction of the whole sample when the pulse starts 
% waveform      - option for the form of the pulse:
%                   - 'sinus' a sine wave,up first (DEFAULT)
%                   - 'box' a boxcar function
%                   - 'gauss' a gaussian pulse
%   NOT YET         - 'raisedcos' raised cosine
% rolloff      - for raised cosine and gauss function, need rolloff factor
%
%OUTPUT:
% xx            - a column vector with the impulse time series and zeros

if nargin < 5
    waveform = 'sinus';
end
if nargin < 6
    startfraction = 0.5*(1-(impulselength/samplength)); % default is halfway
end
if nargin < 7
    rolloff = 0.5;
end

nsamps=samplength/dt;
nsampsw=impulselength/dt;
xx = zeros(nsamps,1); % background trace of zeros

if strcmp(waveform,'sinus')==1
    yy = sin([0:dt:impulselength]*2*pi/impulselength);
elseif strcmp(waveform,'box')==1
    yy = ones(nsampsw+1,1);
elseif strcmp(waveform,'gauss')==1
    tt = dt*(0:nsampsw)';
    tt = tt - (0.5*impulselength);
    sigma = impulselength./8;
    yy = exp((-1/(2*sigma^2))*tt.^2);
else
    error('That pulse shape not supported yet')
end
startind = round(nsamps*startfraction);
xx(startind : startind + nsampsw) = amp * yy; 

end

