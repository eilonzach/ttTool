% produces structure "S" which contains the paths and amplitude information
% for the stations (can be any stations, but I assume you'll be doing this
% for the US array)
% 
% produces "slats", "slons", "samps", which have amplitude at each station
%% prelims
pathfile   = strcat(evdir,'path_stat_',Phase);
ofile_stat = sprintf('%s%samplitudes_stat',evdir,Phase);
if justplot==0
fprintf('\nComputing station amps for phase %s\n',Phase);
%% STA details 
fid = fopen(pathfile,'r');
nsta = cell2mat(textscan(fid,'#NSTA: %f\n',1));
% stas = [1:nsta]'; % will need to change this to list of sta names
slons = zeros(nsta,1);
slats = zeros(nsta,1);
samps = zeros(nsta,1);
%% Set up structure for different stations and phases
S = struct([]);
tic
for is = 1:nsta

if mod(is-1,round(nsta/20)) < 1; fprintf('%g%% done...\n',(100/20)*round(20*is/nsta)); end % progress report
%% read in sta info
sedata=cell2mat(textscan(fid,'#PATH: %f %f %f %f %f\n',1));
if is==1
    elat = sedata(1);
    elon = sedata(2);
    edep = sedata(3);
end
S(is).sta = is; %will need to change this to char(stas(is)) when we have names
S(is).slat = sedata(4);
S(is).slon = sedata(5);
slats(is) = S(is).slat;
slons(is) = S(is).slon;
% read rest of lines
model=textscan(fid,'#MODEL:%s\n',1);
ifani=cell2mat(textscan(fid,'#IFANI:%f\n',1));
if ifani(1)==0, Vph = Vpv; Vsh = Vsv; end % weird way of doing isotropic velocities
ifsh=textscan(fid,'#IFSH:%f\n',1);
phase=textscan(fid,'#PHASE:%s\n',1);
if psvsh_all==0
    psvsh = ispsvsh(char(phase{1}));
else
    psvsh = psvsh_all;
end
narr = cell2mat(textscan(fid,'#NARR: %f\n',1));
if narr==0; textscan(fid,'%s',1,'delimiter','\n'); continue, end % skip if no arrivals
%% LOOP OVER ARRIVALS
for ia = 1:narr    
%% read in path details
ttpara=cell2mat(textscan(fid,'#TTPARA: %f %f %f %f %f %f',1,'delimiter','\n'));
rpath = struct([]);
rpath(1).N = ttpara(1);
if rpath.N == 0; continue; end
D=textscan(fid,'%f%f%f%f%f%f%f%f',rpath.N,'delimiter','\n'); 
rpath.i = D{1};
rpath.rayt = D{2}; if psvsh == 3, rpath.rayt = rpath.rayt./3; end %make SH if needed
rpath.ilay = D{3};
rpath.r = D{4};
rpath.gcd = D{5};
rpath.lat = D{6};
rpath.lon = D{7};
rpath.tt = D{8};
%work out ilay for each rpath
[rpath.ilay] = find_ilay(rpath.r,Rb);
%% calc path details
rmin = min(rpath.r); % minimum r - for calculating ray param
dt = rpath.tt(2:end)-rpath.tt(1:end-1); % time for each sections
% indmin = find(rpath.r==rmin);
rayp = r2d(ttpara(3))/1000;%rmin/linterp(R,Vpv,rmin);
Vpv0 = linterp(R,Vpv,Re-edep); Vph0 = linterp(R,Vph,Re-edep); 
Vsv0 = linterp(R,Vsv,Re-edep); Vsh0 = linterp(R,Vsh,Re-edep);
rho0 = linterp(R,rho,Re-edep);
%incidence angle
if psvsh == 1
inc = asind(rayp*Vpv0/(Re-edep)); % as Re-edep is radius of event
elseif psvsh == 2
inc = asind(rayp*Vsv0/(Re-edep)); % as Re-edep is radius of event
elseif psvsh == 3
inc = asind(rayp*Vsh0/(Re-edep)); % as Re-edep is radius of event
end
faz = azimuth(elat,elon,S(is).slat,S(is).slon);
baz = azimuth(S(is).slat,S(is).slon,elat,elon);
%% calc sections
sect = struct([]);
sect(1).rmid = (rpath.r(1:end-1)+rpath.r(2:end))/2;
sect.ds = sqrt(rpath.r(1:end-1).^2 + rpath.r(2:end).^2 - 2*rpath.r(1:end-1).*rpath.r(2:end).*cosd(rpath.gcd(2:end)-rpath.gcd(1:end-1))); % section length, cosine rule
sect.amid = linterp(R,Vpv,sect.rmid); % p-vel
sect.Bmid = linterp(R,Vsv,sect.rmid); % s-vel
sect.Lmid = (4/3)*(sect.Bmid./sect.amid).^2; % define L = (4/3) (B^2/a^2)
isp = find(abs(rpath.rayt(2:end))==2); % find p-wave sections
iss = find(abs(rpath.rayt(2:end))~=2); % find s-wave sections
% find Q-values in sections (Qa or Qb where appropriate)
sect.Q = zeros(size(sect.rmid));
sect.Q(iss) = linterp(R,Qm,sect.rmid(iss)); %Q_B = Qm
sect.Q(isp) = ( sect.Lmid(isp)./linterp(R,Qm,sect.rmid(isp)) + (1-sect.Lmid(isp))./linterp(R,Qm,sect.rmid(isp)) ).^-1; %Q_a = 1/[L/Qm +(1-L)Qk]
sect.Aq = exp(-0.5*w0*dt./sect.Q);
%% calc piercing
% NB piercing pt assignment will fail if one section passes through more
% than one boundary...
% assign as 999 so can delete later if still 999.
pr = 999*ones(rpath.N,1);
ptt = 999*ones(rpath.N,1);
plon = 999*ones(rpath.N,1);
plat = 999*ones(rpath.N,1);
pgcd = 999*ones(rpath.N,1);
pind = 999*ones(rpath.N,1);
for ipp = 2:Nlay
[ind,x] = crossing(rpath.r,rpath.r,Rb(ipp));	pr(ind) = x; pind(ind) = ind;
[~,x] = crossing(rpath.r,rpath.tt,Rb(ipp));     ptt(ind) = x;
[~,x] = crossing(rpath.r,rpath.lat,Rb(ipp));	plat(ind) = x;
[~,x] = crossing(rpath.r,rpath.lon,Rb(ipp));	plon(ind) = x;
[~,x] = crossing(rpath.r,rpath.gcd,Rb(ipp));    pgcd(ind) = x;
end
% add ends
pr([1;end])   = rpath.r([1;end]);     
ptt([1;end])  = rpath.tt([1;end]);
plat([1;end]) = rpath.lat([1;end]);
plon([1;end]) = rpath.lon([1;end]);
pgcd([1;end]) = rpath.gcd([1;end]);
pind([1;end]) = [1;rpath.N];     
pierce = struct([]);
pierce(1).r = pr(pr~=999);
pierce.tt   = ptt(ptt~=999);
pierce.lon	= plon(plon~=999);
pierce.lat	= plat(plat~=999);
pierce.gcd	= pgcd(pgcd~=999);
pierce.ind  = pind(pgcd~=999);
%% calc bouncing & turning & diffraction
[indbt,rbt] = crossing(gradient(rpath.r,rpath.gcd),rpath.r); % indices of both bouncing and turning pts
indt = indbt(rpath.r(indbt)-rpath.r(indbt-ones(size(indbt)))<0); % indices of turning pts
indd = indbt(rpath.r(indbt)-rpath.r(indbt-ones(size(indbt)))==0); % indices of diffracted pts
indb = indbt(rpath.r(indbt)-rpath.r(indbt-ones(size(indbt)))>0); % indices of bouncing pts
% assign topside reflections from turn to bounce
if isempty(intersect(Rb,rpath.r(indt)))~=1
indtr = find(rpath.r(indt) == intersect(Rb,rpath.r(indt)));
indb = horzcat(indb,indt(indtr)); sort(indb);
indt(indtr) = [];
end
bounce = struct([]);
bounce(1).r = rpath.r(indb);
bounce.tt  = rpath.tt(indb);
bounce.lon = rpath.lon(indb);
bounce.lat = rpath.lat(indb);
bounce.gcd = rpath.gcd(indb);
diffr = struct([]); 
diffr(1).r = rpath.r(indd);
diffr.tt  = rpath.tt(indd);
diffr.lon = rpath.lon(indd);
diffr.lat = rpath.lat(indd);
diffr.gcd = rpath.gcd(indd);
turn = struct([]); 
turn(1).r = rpath.r(indt);
turn.tt  = rpath.tt(indt);
turn.lon = rpath.lon(indt);
turn.lat = rpath.lat(indt);
turn.gcd = rpath.gcd(indt);
%% Reflection Transmission coefficients
% Transmission
TC = zeros(length(pierce.r)/2 - 1,1);
for itt = 1:length(TC)
    indtr = find(R == pierce.r(2*itt));
    a = Vpv(indtr);
    B = Vsv(indtr); B(B==0)=1e-10;
    d = rho(indtr);
    TC(itt) = conversion_coefficient(rpath.rayt(pierce.ind(2*itt)),rpath.rayt(pierce.ind(2*itt + 1)),...
            rayp/pierce.r(2*itt),a(2),B(2),d(2),a(1),B(1),d(1));
end
% Reflection
RC = zeros(length(bounce.r),1);
for irr = 1:length(RC)
    indr = find(R == bounce.r(irr));
    a = Vpv(indr); if length(indr)==1, a = [a,0]; end
    B = Vsv(indr); if length(indr)==1, B = [B,0]; end
    d = rho(indr); if length(indr)==1, d = [d,0]; end
    RC(irr) = conversion_coefficient(rpath.rayt(indb(irr)),rpath.rayt(indb(irr)+1),...
            rayp/bounce.r(irr),a(end),B(end),d(end),a(1),B(1),d(1));
end
%% Diffraction effect...
diffr.ds = sqrt(diffr.r(1:end-1).^2 + diffr.r(2:end).^2 - 2*diffr.r(1:end-1).*diffr.r(2:end).*cosd(diffr.gcd(2:end)-diffr.gcd(1:end-1))); % section length, cosine rule
AD = exp(-sum(100*rayp*diffr.ds/Re));
%% Free surface correction
W = freesurf(rayp/Re,Vpv(end),Vsv(end));
% think this is a fudge! 
if psvsh == 1
    W = abs(W(1,1));
elseif psvsh == 2
    W = abs(W(2,2));
elseif psvsh == 3
    W = 2;
end
%% Amplitudes
ATRW = abs(prod(TC)*prod(RC))*W;
AQ = prod(sect.Aq);
% displacements
s = sum(sect.ds); %total rpathlength (km)
As = 1/(s*1000); % geometric spreading
% radiation pattern
[ Rpat ] = radpcalc(str,dip,rak,faz,inc);
if psvsh == 1
    A0 = 1./(4*pi*rho0*Vpv0^3);
elseif psvsh ==2
    A0 = 1./(4*pi*rho0*Vsv0^3);
elseif psvsh ==3
    A0 = 1./(4*pi*rho0*Vsh0^3);
end
A = A0*ATRW*As*AQ*Rpat(psvsh)*AD; % final magnitude
%% results structures for this phase
phs = struct([]);
phs(1).phase = char(phase{1});
phs.rayp = rayp;
phs.rmin = rmin;
phs.tt = rpath.tt(end);
phs.gcd = rpath.gcd(end);
phs.inc = inc;
phs.faz = faz;
phs.baz = baz;
phs.psvsh = psvsh;
phs.AQ = AQ;
phs.AS = As;
phs.A = A;
phs.rpath = rpath;
phs.sect = sect;
phs.pierce = pierce;
phs.turn = turn;
phs.bounce = bounce;
phs.diffr = diffr;
S(is).phase = char(phase{1});
S(is).rays(ia) = phs;
end % arrivals loop

samps(is) = S(is).rays(1).A;
% %% fudge caustic...
% if length(S(is).rays) > 1
%     if S(is).rays(2).tt - S(is).rays(1).tt < 15
%         samps(is) = samps(is) + S(is).rays(2).A;
%     end
% end
end % stas loop
fprintf('100%% done... only took %.3f seconds\n',toc)
fclose(fid);
end % if just plot
