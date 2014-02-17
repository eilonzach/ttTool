function [ mrlR,mu,Var,sdev ] = mrl( theta,r,option )
% [ mrlR,mu,Var,sdev ] = mrl( theta,length,option )
% calculates the mean resultant length of a set of angles in RADIANS
% the option is for 'axial' (i.e. 181=1) or 'direct' data
% NB can omit r
% OUTPUTS
% mrlR  - the mean resultant length of the average vector
% mu    - the direction of the average vector
% Var   - the sample circular variance
% sdev  - the sample circular standard deviation
%
% EDIT 26/10/12
% can now take matrices theta, r - computes mrlR etc over columns
% i.e. 10x4 theta matrix will give four 1x4 row matrices of mu, mrlR, etc.
%

if nargin < 2
    r = ones(size(theta));
end

if nargin == 2 % if you put just angles and then option for axial or direct
    if ischar(r)==1
        option = r;
        r = ones(size(theta));
    end
end

if nargin == 3
    if strcmp(option,'axial')==1
        theta=mod(2*theta,2*pi);
        theta=mp2pp(theta);
    elseif strcmp(option,'direct')==1
        theta=mp2pp(theta);
    end
end   

n = size(theta,1);
m = size(theta,2);
mrlR = zeros(1,m);
mu   = zeros(1,m);
Var  = zeros(1,m);
sdev = zeros(1,m);

for im = 1:m
% fprintf('number of samples, n, is %i \n',n);
C=sum(r(:,im).*cos(theta(:,im)));
S=sum(r(:,im).*sin(theta(:,im)));
R=sqrt(C^2 + S^2);
mrlR(im)=R/n;
if C==0; error('angles sum to nothing'); end
if S==0; mu=0; end
% fprintf('sample mean resultant length, mrlR, is %.3f \n',mrlR);
if S>0 && C>0
    mu(im)=atan(S/C); % in radians
elseif S<0 && C>0
    mu(im)=atan(S/C) + 2*pi; % in radians
elseif C<0
    mu(im)=atan(S/C)+pi; % in radians
end
end

if nargin==3
    if strcmp(option,'axial')==1
        mu=0.5*mu;
    end
end   

% Var = 1-mrlR;
Var = 1-(mrlR./mean(r));
sdev=sqrt(-2.*log(1-Var));

% mud=r2d(mu); % in degrees
end

