% Makes relevant plots for the travel time and amplitude plots


%% plotting prereqs.
tic
figure(89), clf; 

latlim = [22 50];
lonlim = [-130 -70];
lonnodes = [lonlim(1):0.5:lonlim(2)];
latnodes = [latlim(1):0.5:latlim(2)];
cc = 10;

[ampgrid,longrid,latgrid] = gridfit(glons,glats,gamps,lonnodes,latnodes);
ampf = max(max(abs(ampgrid)));

[gcarcmin,~]   = min(distance(elat,elon,slats,slons));
[gcarcmax,ind] = max(distance(elat,elon,slats,slons));
scrcs = cc*[ceil(gcarcmin/cc):floor(gcarcmax/cc)]';
faz = azimuth(elat,elon,slats(ind),slons(ind));
[midlat,midlon] = reckon(elat,elon,gcarcmax/2,faz);
[lablat,lablon] = reckon(elat,elon,scrcs,faz);

psvsh = ispsvsh(Phase);

%% 2) Radiation pattern
radpplot(str,dip,rak,psvsh,89,[2 3 3]);
h2 = gca; p2 = get(h2,'pos');

hold on
for is = 1:nsta
if isempty(S(is).rays), continue, end
[ x,y ] = Lowhemis_conv(S(is).rays(1).faz,90-S(is).rays(1).inc,4);
plot(x,y,'.g')
end
hold off
freezeColors

%% 3) Path map
figure(89); 
h3 = subplot(2,3,4); p3 = get(h3,'pos'); % get position of axes
set(h3,'pos',[p3(1)-0.05 p3(2)-0.04 p3(3)+0.08 p3(4)+0.08])

m_proj('UTM','longitudes',lonlim,'latitudes',latlim);
% m_gshhs_l('patch',[.9 .9 .9],'edgecolor','k');
m_coord('geographic');
m_grid;
hold on
for is = 1:nsta
    if isempty(S(is).rays), continue, end
    m_plot(S(is).rays(1).rpath.lon,S(is).rays(1).rpath.lat,'k')
    m_plot(S(is).rays(1).pierce.lon,S(is).rays(1).pierce.lat,'.k')
end
 m_range_ring(elon,elat,Re*d2r([0:cc:180]),'color','b','linewi',2);
 h = m_text(lablon,lablat,cellstr(num2str(scrcs)),'clip','on');
 set(h,'BackgroundColor','w','Margin',1,'color','b','FontWeight','bold','HorizontalAlignment','center')
%  h = m_text(lablon,lablat,cellstr(num2str([0:cc:180]')),'clip','on');
%  set(h,'color','w','HorizontalAlignment','center')

 % smallcircs = zeros(181,36,2); % each column is a different distance, first layer is lats second layer is lons
% for igc = 1:36
%     [la,lo] = reckon(elat,elon,5*igc,[0:2:360]');
%     [x,y] = m_ll2xy(lo,la,'clip','on');
%     [smallcircs(:,igc,2),smallcircs(:,igc,1)] = m_xy2ll(x,y);
% end
% [x,y] = m_ll2xy(smallcircs(:,:,2),smallcircs(:,:,1),'clip','on');
% [la,lo] = m_xy2ll(x,y)
% m_plot(smallcircs(:,:,2),smallcircs(:,:,1),'k')
hold off

% ht3 = title('Paths','FontSize',14); pt3 = get(ht3,'pos');

%% 4) plot rayrpath in side view
figure(89); 
h4 = subplot(2,3,5); p4 = get(gca,'pos');
set(h4,'pos',[p4(1)+0.005, p4(2)-0.01 p4(3)+0.01 p4(4)+0.01])

th = linspace(0,2*pi,721)';
% polar2(th*ones(size(Rb))',ones(size(th))*Rb',[0 Re])
hold on
for il = 1:Nlay
polar2(th,ones(size(th))*Rb(il),[0 Re]);
end
for is = 1:nsta
    if isempty(S(is).rays), continue, end
    if is ~= 20*round(is/20), continue; end % only do every 20th
polar2(d2r(S(is).rays(1).rpath.gcd),S(is).rays(1).rpath.r,[0 Re],'r');
polar2(d2r(S(is).rays(1).pierce.gcd),S(is).rays(1).pierce.r,[0 Re],'r');
polar2(d2r(S(is).rays(1).turn.gcd),S(is).rays(1).turn.r,[0 Re],'.r');
polar(d2r(S(is).rays(1).bounce.gcd),S(is).rays(1).bounce.r,'.g');
end
hold off

ht4 = title('Raypaths side slice','FontSize',14); pt4 = get(ht4,'pos');
set(ht4,'pos',[-7.7e3 0 0]);
view(90,-90)
axis equal; box off; axis off

%% 5) Full Path map
figure(89); 
h5 = subplot(2,3,6); p5 = get(h5,'position');

m_proj('ortho','lat',midlat','long',midlon');
m_coast('patch',[0.9 0.9 0.9]);
m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',[-80:40:80]);
hold on
for is = 1:nsta
if isempty(S(is).rays), continue, end
m_plot(S(is).rays(1).rpath.lon,S(is).rays(1).rpath.lat,'k')
m_plot(S(is).rays(1).turn.lon,S(is).rays(1).turn.lat,'.g')
m_plot(S(is).rays(1).pierce.lon(S(is).rays(1).pierce.r==Rb(3)),...
    S(is).rays(1).pierce.lat(S(is).rays(1).pierce.r==Rb(3)),'.r')
m_plot(S(is).rays(1).bounce.lon,S(is).rays(1).bounce.lat,'.m')
end
m_range_ring(elon,elat,Re*d2r([0:cc:180]),'color','b');
hold off

ht5 = title('Great circle paths','visible','on','FontSize',14); pt5 = get(ht5,'pos'); 
set(ht5,'pos',[0 -1.2 0]);

%% Amplitude map
figure(89),
h1 = subplot(2,3,1:2); p1 = get(h1,'pos'); % get position of axes
set(h1,'pos',[p1(1)-0.1 p1(2)-0.2 p1(3)+0.1 p1(4)+0.3]) % move the axes slightly

% define projection
m_proj('UTM','longitudes',lonlim,'latitudes',latlim);
% m_gshhs_l('patch',[.9 .9 .9],'edgecolor','k');
m_coord('geographic');
m_grid;
hold on

m_pcolor(longrid,latgrid,ampgrid);  % plot amplitudes

shading interp
m_gshhs_l('LineWidth',2,'Color','k');
m_plot(slons,slats,'m^','MarkerFaceColor',[.9 .9 .9],'MarkerSize',4)  % plot stas
hold off
colormap('default')
colormap(flipud(colormap))
caxis(ampf*[-1 1])
hc1 = cbfreeze(colorbar); pc1 = get(hc1,'pos');
set(hc1,'pos',[pc1(1) pc1(2)+0.1 pc1(3) pc1(4)-0.15]);
freezeColors

ht1 = title('Predicted amplitude (relative)','FontSize',14); pt1 = get(ht1,'pos'); 
set(ht1,'pos',[pt1(1)+0.01 pt1(2) pt1(3)]);

fprintf('Plotting took %.3f seconds\n',toc);

%% Amplitude-o-gram
figure(90), clf, 
if ~isempty(refsta), subplot(3,1,1); end
% h = subplot(2,3,1:2); p = get(h,'pos'); % get position of axes
% set(h,'pos',[p(1)-0.1 p(2)-0.2 p(3)+0.1 p(4)+0.3]) % move the axes slightly
title('Predicted amplitude (relative)','FontSize',14);
samplength = 4000;
impulselength = 20;
dt = 0.5;
xlimits = [000,2800];

xx = zeros(samplength./dt,4*nphase);
tt = [0:dt:samplength-dt];
for ip = 1:nphase
    for ia = 1:length(T(ip).rays)
    xx(:,ia + 4*(ip-1)) = synthtrace(samplength,impulselength,T(ip).rays(ia).A,dt,'gauss',T(ip).rays(ia).tt/samplength,dt);
    end
end
for ip = 1:nphase
    if isodd(T(ip).phase(1)); f = -1; else f = 1; end
    for ia = 1:length(T(ip).rays)
        a = max(max(xx));
        y = 1.1*a + f*0.06*a*length(T(ip).phase);
    text(T(ip).rays(ia).tt+impulselength/2,y,phases(ip),'HorizontalAlignment','center');
    end
end

hold on
    plot(tt,abs(sum(xx,2)),'b','LineWidth',2);
    plot(tt,abs(xx),'r');
hold off
xlim(xlimits)
set(gca,'XTick',[],'YTick',[],'ycolor',get(gcf,'color'));

% refsta
if ~isempty(refsta)
ymax = max(double([datr(:);datt(:);datz(:)]));
subplot(3,1,2:3)
plot(sampletimes,datt+1.7*ymax,'g','LineWidth',1.5);
hold on
plot(sampletimes,datr,'r','LineWidth',1.5);          
plot(sampletimes,datz-1.7*ymax,'b','LineWidth',1.5); 
t = text((xlimits(1)-80)*ones(3,1),1.7*ymax.*[1 0 -1]',{'TRANSVERSE';'RADIAL';'VERTICAL'});
set(t, 'rotation', 90,'HorizontalAlignment','center','FontSize',14)
hold off
xlim(xlimits)
xlabel('Seconds since event','Fontsize',14)
end

%% plot rayrpath in side view
figure(91); clf
title('Raypaths side slice','FontSize',14)

th = linspace(0,2*pi,721)';
% polar2(th*ones(size(Rb))',ones(size(th))*Rb',[0 Re])
hold on
for il = 1:Nlay
polar2(th,ones(size(th))*Rb(il),[0 Re],'k');
end
view(90,-90)
for ip = 1:nphase
    for ia = 1:length(T(ip).rays)
        rpath = T(ip).rays(ia).rpath;
        isp = find(abs(rpath.rayt(1:end))==2); % find s-wave sections
        iss = find(abs(rpath.rayt(1:end))~=2); % find p-wave sections
        ppath = rpath;  ppath.gcd(iss) = NaN;  ppath.r(iss) = NaN;
        if T(ip).rays(ia).psvsh == 1
            h1 = polar2(d2r(rpath.gcd),rpath.r,[0 Re],'r'); % plot all as S
            h2 = polar2(d2r(ppath.gcd),ppath.r,[0 Re],'b'); %plot over any P
        else
            h1 = polar2(2*pi-d2r(rpath.gcd),rpath.r,[0 Re],'r'); % plot all as S
            h2 = polar2(2*pi-d2r(ppath.gcd),ppath.r,[0 Re],'b'); %plot over any P
        end
        set([h1,h2],'Linewidth',1.5)
    end
end
for ip = 1:nphase
    for ia = 1:length(T(ip).rays)
        rpath = T(ip).rays(ia).rpath;
        i = random('unid',rpath.N);
        x = rpath.r(i)*sind(rpath.gcd(i));
        y = rpath.r(i)*cosd(rpath.gcd(i));
        if T(ip).rays(ia).psvsh == 1
            text(y,x,T(ip).phase,'FontSize',12,'FontWeight','bold');
        else
            text(y,-x,T(ip).phase,'FontSize',12,'FontWeight','bold');
        end
    end
end
% %% TEMP ADDITION
% h3 = polar2(d2r(T1Hz(1).rays(1).rpath.gcd),T1Hz(1).rays(1).rpath.r,[0 Re],'r');
% set(h3,'Linewidth',1.5)
% %% END TEMP ADDITION
hold off
h = gca;
position = get(h,'Position');
outposition = get(h,'Outerposition');
set(h,'Position', [0.05 0.05 0.9 0.9]);
axis equal; box off; axis off

%% save to pdf
if savefig == 1
h = figure(89);
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(h,'-zbuffer','-dpdf','-r200',strcat(evdir,sprintf('%samp_map',Phase)))

% save2pdf(89,sprintf('%samplitudes_stas',char(phase{1})),evdir);
end