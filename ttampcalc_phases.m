% produces structure "T" which contains the paths and amplitude information
% for each different phase arriving at the reference station
%% prelims
pathfile  = strcat(evdir,'path__1sta');
ofile     = sprintf('%samplitude_o_gram',evdir);
if justplot==0
fprintf('\nComputing 1sta amps\n');
%% STA details 
fid = fopen(pathfile,'r');
nphase = cell2mat(textscan(fid,'#NPHASE: %f\n',1));
%% Set up structure for different stations and phases
T = struct([]);
phases = {};
tic
for ip = 1:nphase
%% read in sta info
sedata=cell2mat(textscan(fid,'#PATH: %f %f %f %f %f\n',1));
if ip==1
    elat = sedata(1);
    elon = sedata(2);
    edep = sedata(3);
    T(1).slat = sedata(4);
    T(1).slon = sedata(5);
end
% read rest of lines
model=textscan(fid,'#MODEL:%s\n',1);
ifani=cell2mat(textscan(fid,'#IFANI:%f\n',1));
if ifani(1)==0, Vph = Vpv; Vsh = Vsv; end % weird way of doing isotropic velocities
ifsh=textscan(fid,'#IFSH:%f\n',1);
phase=textscan(fid,'#PHASE:%s\n',1);
phases = cat(1,phases,phase{1}); phase = char(phase{1});
psvsh = ispsvsh(phase);
narr = cell2mat(textscan(fid,'#NARR: %f\n',1));
T(ip).sta = 'station';
T(ip).phase = phase;
if narr==0; textscan(fid,'%s',1,'delimiter','\n'); continue, end % skip if no arrivals
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
faz = azimuth(elat,elon,T(1).slat,T(1).slon);
baz = azimuth(T(1).slat,T(1).slon,elat,elon);
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
% A = A0*As*AQ*Rpat(psvsh)*W; % final magnitude
%% results structures for this phase
phs = struct([]);
phs(1).phase = phase;
phs.rayp = rayp;
phs.rmin = rmin;
phs.tt = rpath.tt(end);
phs.gcd = rpath.gcd(end);
phs.inc = inc;
phs.faz = faz;
phs.baz = baz;
phs.psvsh = psvsh;
phs.Rpat = Rpat;
phs.ATR = ATRW;
phs.AQ = AQ;
phs.AS = As;
phs.A = A;
phs.rpath = rpath;
phs.sect = sect;
phs.pierce = pierce;
phs.turn = turn;
phs.bounce = bounce;
phs.diffr = diffr;
T(ip).rays(ia) = phs;

% if strcmp(phase{1},'P'); return; end
end % arrivals loop

% samps(ip) = T(ip).rays(1).A;
% %% fudge caustic...
% if length(T(ip).rays) > 1
%     if T(ip).rays(2).tt - T(ip).rays(1).tt < 15
%         samps(ip) = samps(ip) + T(ip).rays(2).A;
%     end
% end
end % phases loop
fprintf('100%% done... only took %.3f seconds\n',toc)
fclose(fid);

%% Fetch data from IRIS
if ~isempty(refsta)
%   args: Net, Sta, Loc, Cha, Starttime, Endtime [,quality][,includePZ][,verbosity]
t0str = sprintf('%u-%u-%u %u:%u:%u',evtime);
t1str = sprintf('%u-%u-%u %u:%u:%u',evtime + [0 0 0 1 0 0]);
mytrace = irisFetch.Traces('IU',refsta,'*','BH*',t0str,t1str,'includePZ');
tracee = mytrace(1);
tracen = mytrace(2);
tracez = mytrace(3);
sampletimes=[0:tracee.sampleCount-1]./tracee.sampleRate;
%% Rotate data to get R,T
datr = [tracen.data,tracee.data]*[cosd(T(1).rays(1).baz + 180);sind(T(1).rays(1).baz + 180)];
datt = [tracen.data,tracee.data]*[-sind(T(1).rays(1).baz + 180);cosd(T(1).rays(1).baz + 180)];
datz = tracez.data;
%% take gradient and then negative to get into (scaled) displacement
% datr = -gradient(datr,1./tracee.sampleRate);
% datt = -gradient(datt,1./tracee.sampleRate);
% datz = -gradient(datz,1./tracee.sampleRate);
end % get data if refsta given

end % if just plot

%% Print travel times
fprintf('Travel times for station at %.1f degrees:\n',distance(elat,elon,T(1).slat, T(1).slon))
for ip = 1:nphase
    for ia = 1:length(T(ip).rays)
        if ia==1
        fprintf('%s\t%7.2f s\n',T(ip).phase,T(ip).rays(ia).tt)
        else
        fprintf('     \t%7.2f s\n',T(ip).rays(ia).tt)
        end
    end % arr loop
end % phase loop

% %% plotting prereqs.
% tic
% 
% %% plot rayrpath in side view
% figure(89); clf
% title('Raypaths side slice','FontSize',14)
% 
% th = linspace(0,2*pi,721)';
% % polar2(th*ones(size(Rb))',ones(size(th))*Rb',[0 Re])
% hold on
% for il = 1:Nlay
% polar2(th,ones(size(th))*Rb(il),[0 Re],'k');
% end
% view(90,-90)
% for ip = 1:nphase
%     for ia = 1:length(T(ip).rays)
%         rpath = T(ip).rays(ia).rpath;
%         isp = find(abs(rpath.rayt(1:end))==2); % find s-wave sections
%         iss = find(abs(rpath.rayt(1:end))~=2); % find p-wave sections
%         ppath = rpath;  ppath.gcd(iss) = NaN;  ppath.r(iss) = NaN;
%         if T(ip).rays(ia).psvsh == 1
%             h1 = polar2(d2r(rpath.gcd),rpath.r,[0 Re],'r'); % plot all as S
%             h2 = polar2(d2r(ppath.gcd),ppath.r,[0 Re],'b'); %plot over any P
%         else
%             h1 = polar2(2*pi-d2r(rpath.gcd),rpath.r,[0 Re],'r'); % plot all as S
%             h2 = polar2(2*pi-d2r(ppath.gcd),ppath.r,[0 Re],'b'); %plot over any P
%         end
%         set([h1,h2],'Linewidth',1.5)
%     end
% end
% for ip = 1:nphase
%     for ia = 1:length(T(ip).rays)
%         rpath = T(ip).rays(ia).rpath;
%         i = random('unid',rpath.N);
%         x = rpath.r(i)*sind(rpath.gcd(i));
%         y = rpath.r(i)*cosd(rpath.gcd(i));
%         if T(ip).rays(ia).psvsh == 1
%             text(y,x,T(ip).phase,'FontSize',12,'FontWeight','bold');
%         else
%             text(y,-x,T(ip).phase,'FontSize',12,'FontWeight','bold');
%         end
%     end
% end
% % %% TEMP ADDITION
% % h3 = polar2(d2r(T1Hz(1).rays(1).rpath.gcd),T1Hz(1).rays(1).rpath.r,[0 Re],'r');
% % set(h3,'Linewidth',1.5)
% % %% END TEMO ADDITION
% hold off
% h = gca;
% position = get(h,'Position');
% outposition = get(h,'Outerposition');
% set(h,'Position', [0.05 0.05 0.9 0.9]);
% axis equal; box off; axis off
% 
% %% Radiation pattern
% figure(91); clf
% radpplot(str,dip,rak,1,91,[1,3,1]);
% radpplot(str,dip,rak,2,91,[1,3,2]);
% radpplot(str,dip,rak,3,91,[1,3,3]);
% for ip = 1:nphase
%     if isempty(T(ip).rays), continue, end
%     for ia = 1:length(T(ip).rays)
%         r = 4 * (sind(T(ip).rays(ia).inc/2)/sind(45));
%         x = sind(T(ip).rays(ia).faz)*r';
%         y = cosd(T(ip).rays(ia).faz)*r';
%         subplot(1,3,T(ip).rays(ia).psvsh)
%         hold on
%         plot(x,y,'.k')
%         text(x,y+0.5,T(ip).rays(ia).phase,'HorizontalAlignment','center')
%         hold off
%     end
% end
% fprintf('Plotting took %.3f seconds\n',toc);
% freezeColors
% 
% %% Times
% figure(92), clf, 
% % h = subplot(2,3,1:2); p = get(h,'pos'); % get position of axes
% % set(h,'pos',[p(1)-0.1 p(2)-0.2 p(3)+0.1 p(4)+0.3]) % move the axes slightly
% title('Arrival times','FontSize',14);
% xlimits = [700,2500];
% ylimits = [-4e5 4e5];
% ymax = max(double([datr(:);datt(:);datz(:)]));
% hold on
% for ip = 1:nphase
%     if isodd(T(ip).phase(1)); f = -1; else f = 1; end
%     for ia = 1:length(T(ip).rays)
%         line(T(ip).rays(ia).tt*[1,1],2.3*ymax*[-1 1])
%         if ia==1
%         t = text(T(ip).rays(1).tt,      f*(2.3+0.12*length(T(ip).phase))*ymax,phases(ip));
%         set(t,'HorizontalAlignment','center','FontSize',10)
%         end
%     end
% end
% plot(sampletimes,datt+1.7*ymax,'r','LineWidth',1.5);
% plot(sampletimes,datr,'r','LineWidth',1.5);          
% plot(sampletimes,datz-1.7*ymax,'r','LineWidth',1.5); 
% t = text((xlimits(1)-50)*ones(3,1),1.7*ymax.*[1 0 -1]',{'TRANSVERSE';'RADIAL';'VERTICAL'});
% set(t, 'rotation', 90,'HorizontalAlignment','center','FontSize',14)
% hold off
% set(gca,'YTick',[],'ycolor',get(gcf,'color'));
% xlim(xlimits)
% ylim(ylimits)
% xlabel('Seconds since event','Fontsize',14)
% return
% %% Amplitude-o-gram
% figure(90), clf, 
% subplot(3,1,1)
% % h = subplot(2,3,1:2); p = get(h,'pos'); % get position of axes
% % set(h,'pos',[p(1)-0.1 p(2)-0.2 p(3)+0.1 p(4)+0.3]) % move the axes slightly
% title('Predicted amplitude (relative)','FontSize',14);
% samplength = 4000;
% impulselength = 20;
% dt = 0.5;
% xlimits = [700,2800];
% ymax = max(double([datr(:);datt(:);datz(:)]));
% 
% xx = zeros(samplength./dt,4*nphase);
% tt = [0:dt:samplength-dt];
% for ip = 1:nphase
%     for ia = 1:length(T(ip).rays)
%     xx(:,ia + 4*(ip-1)) = synthtrace(samplength,impulselength,T(ip).rays(ia).A,dt,'gauss',T(ip).rays(ia).tt/samplength,dt);
%     end
% end
% for ip = 1:nphase
%     if isodd(T(ip).phase(1)); f = -1; else f = 1; end
%     for ia = 1:length(T(ip).rays)
%         a = max(max(xx));
%         y = 1.1*a + f*0.06*a*length(T(ip).phase);
%     text(T(ip).rays(ia).tt+impulselength/2,y,phases(ip),'HorizontalAlignment','center');
%     end
% end
% 
% hold on
%     plot(tt,abs(sum(xx,2)),'b','LineWidth',2);
%     plot(tt,abs(xx),'r');
% hold off
% xlim(xlimits)
% set(gca,'XTick',[],'YTick',[],'ycolor',get(gcf,'color'));
% 
% subplot(3,1,2:3)
% plot(sampletimes,datt+1.7*ymax,'g','LineWidth',1.5);
% hold on
% plot(sampletimes,datr,'r','LineWidth',1.5);          
% plot(sampletimes,datz-1.7*ymax,'b','LineWidth',1.5); 
% t = text((xlimits(1)-80)*ones(3,1),1.7*ymax.*[1 0 -1]',{'TRANSVERSE';'RADIAL';'VERTICAL'});
% set(t, 'rotation', 90,'HorizontalAlignment','center','FontSize',14)
% hold off
% xlim(xlimits)
% xlabel('Seconds since event','Fontsize',14)
% 
% %% save to pdf
% if savefig == 1
% h = figure(90);
% set(h,'PaperOrientation','landscape');
% set(h,'PaperUnits','normalized');
% set(h,'PaperPosition', [0 0 1 1]);
% print(h,'-zbuffer','-dpdf','-r200',ofile)
% h = figure(91);
% set(h,'PaperOrientation','landscape');
% set(h,'PaperUnits','normalized');
% set(h,'PaperPosition', [0 0 1 1]);
% print(h,'-zbuffer','-dpdf','-r200',strcat(evdir,'rpat'))
% h = figure(89);
% set(h,'PaperOrientation','landscape');
% set(h,'PaperUnits','normalized');
% set(h,'PaperPosition', [0 0 1 1]);
% print(h,'-zbuffer','-dpdf','-r200',strcat(evdir,'raypaths'))
% % save2pdf(89,sprintf('%samplitudes_stas',char(phase{1})),evdir);
% end