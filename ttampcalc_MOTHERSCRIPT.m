% Mother script for all the travel time and amplitude calculation tools
%% prelims
justplot = 0; % option to skip calculations and just re-plot
savefig = 1; % option so save output figure as pdf
psvsh_all = 0; % can choose psvsh - 0 sets to default, using ispsvsh 
w0 = 2*pi/1; % angular freq.
refsta = ''; % name of ref. station - to get s-grams - or ''

%% files - make sure all these paths are correct
modelfile = '/Users/Zach/Documents/MATLAB/TTamp_tool/NOH2O';
evdir     = '/Users/Zach/Documents/MATLAB/TTamp_tool/201401200252/'; % remember a final forward slash!
CMTfile = strcat(evdir,'CMTSOLUTION');
phasesfile =strcat(evdir,'phases');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% NO NEED TO ALTER BELOW HERE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Earthquake details 
eq = CMTSOLUTIONread(strcat(evdir,'CMTSOLUTION'));
[M0, str, dip, rak] = momenttensor_invert(eq.M);
str = str(1);
dip = dip(1);
rak = rak(1);
R = [0,-1,0;0,0,1;-1,0,0];  % rotation matrix from Up,S,E to N,E,Down
M = R*eq.M*R'; 
%% read in model data
fid = fopen(modelfile);
Mf=textscan(fid,'%f%f%f%f%f%f%f%f%f','headerlines',3); 
fclose(fid);
R = Mf{1}/1000; % in km
Re = max(R); % Earth radius
Rb = [0;R(diff(R)==0);Re]; % Radii of boundaries, inc centre and edge
Nlay = length(Rb)-1;
rho = Mf{2};
Vpv = Mf{3};
Vsv = Mf{4};
Qk  = Mf{5};
Qm  = Mf{6}; Qm(Qm==0)=99999999;
Vph = Mf{7};
Vsh = Mf{8};
%% Phases
fid = fopen(phasesfile,'r');
Phases = textscan(fid,'%s'); 
Phases = Phases{1};
fclose(fid);
fprintf('Phases:\n')
for ip = 1:length(Phases)
    fprintf('%.0f: %s\n',ip,char(Phases(ip)))
end
ip = input('Choose phase (number):  ');
Phase = char(Phases(ip));
%% TTampcalc_grid
run('/Users/Zach/Documents/MATLAB/TTamp_tool/ttampcalc_grid.m')
%% TTampcalc_stations
run('/Users/Zach/Documents/MATLAB/TTamp_tool/ttampcalc_stations.m')
%% TTampcalc_phases 
run('/Users/Zach/Documents/MATLAB/TTamp_tool/ttampcalc_phases.m')
%% TTampcalc_plotting
run('/Users/Zach/Documents/MATLAB/TTamp_tool/ttampcalc_plotting.m')

