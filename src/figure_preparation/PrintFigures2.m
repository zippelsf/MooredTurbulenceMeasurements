% Script to generate figures for 2nd round of paper
% -Much of this code was developed in "/Sandbox" before moving here
% -Saves figurese as .pdf as a default, uses shell script to auto-move to my
%   overleaf folder
% - Last modified: SFZ 5/21

clear
close all
clc

%% Figure 2

figure(2),clf
SyntheticSampleTest();
SPURS_figsave(gcf)

%% Figure 3 [formerly fig. 6a] SNR/Noise Floor plot (current fig.3)
clear
close all
addpath('../processing/')
burst = load('../../Data/SPURS-2/interim/NortekFiles/8116/811604_burst0510.mat');
burst_1MHz = load('../../Data/SPURS-1/interim/NortekFiles/9134/913411_burst2879.mat');

figure(3)
clf
set(gcf,'Position',[106   341   837   455])
% set(gcf,'Position',[173   385   864   390])

e = [-9,-8,-7,-6];
for ei = 1:4
%make theory spectrum, including dissipative roll-off
k_grid = 0.5:0.1:1000;
epsilon = 10^e(ei);
nu = 10^-6;
eta = (nu.^3 / epsilon)^(1/4);
beta = 5.2;
c_eta = 0.4;
f_eta = exp(-beta*( ((k_grid*eta).^4 + c_eta^4).^(1/4) -c_eta) );
c1 = 1.5;

%Make one-sided via numeric integration
E_pope_omni = c1*epsilon^(2/3).*k_grid.^(-5/3) .*f_eta;
for ki = 1: (length(k_grid)-1)
    inds = ki:length(k_grid);
    E_pope(ki) = trapz(k_grid(inds), E_pope_omni(inds)./k_grid(inds).*(1 - k_grid(ki)^2./k_grid(inds).^2) );
end
E_pope(length(k_grid)) = NaN;

[k_nasmyth, P_nasmyth]=univ_turb_spec(epsilon); %in rad/m or 1/m??


p1=loglog(k_grid/2/pi,E_pope*2*pi,'-k','linewidth',2 );
hold on
%p2=plot(k_nasmyth, P_nasmyth,'-r','linewidth',2); Changed 4/27/21

end


% calc noise floors, and variance ceilings for 1Mhz and 2 MHz
M = 10;
R = 95;
Vr = burst.Vr;
floor_2MHz = (Vr/pi)^2 * ((R/100).^(-2) - 1)/(2*M)./burst.k(end);
floor2_2MHz = (Vr/pi)^2 * ((60/100).^(-2) - 1)/(2*M)./burst.k(end);

floor_2MHz_zedel = -2*log(60/100)*(Vr/pi)^2./burst.k(end)./(M);

[k_nasmyth, P_nasmyth]=univ_turb_spec(4*10^-4); %in rad/m or 1/m??
[~,ki] = min( abs(burst.k(2)/2/pi - k_nasmyth) );
ceil_2MHz = P_nasmyth(ki)/2/pi;
box_2MHz =  [burst.k(2)/2/pi, burst.k(end)/2/pi ; floor_2MHz*2*pi, ceil_2MHz*2*pi];

clrs = get( groot, 'defaultAxescolororder'); %COLORS
b1=plot(box_2MHz(1,:),box_2MHz(2,[1,1]),'--','linewidth',2,'color',clrs(4,:));
plot(box_2MHz(1,:),floor2_2MHz*[1,1].*2.*pi,'--','linewidth',2,'color',clrs(4,:));
%plot(box_2MHz(1,:),floor_2MHz_zedel*[1,1].*2.*pi,'--','linewidth',2,'color',clrs(4,:));
plot(box_2MHz(1,:),box_2MHz(2,[2,2]),'--','linewidth',2,'color',clrs(4,:))
plot(box_2MHz(1,[1,1]),box_2MHz(2,[1,2]),'--','linewidth',2,'color',clrs(4,:))
plot(box_2MHz(1,[2,2]),box_2MHz(2,[1,2]),'--','linewidth',2,'color',clrs(4,:))
patch2mhz=patch([box_2MHz(1,:),fliplr(box_2MHz(1,:))],...
    [box_2MHz(2,1),box_2MHz(2,1),floor2_2MHz*2*pi,floor2_2MHz*2*pi],...
    clrs(4,:),'FaceAlpha',0.2);

M = 13;
Vr = 0.085;
floor_1MHz = (Vr/pi)^2 * ((R/100).^(-2) - 1)/(2*M)./burst_1MHz.k(end); %exponent changed to R^-2!! Note confusing notation diffs with Zedel and Scherebina papers
floor2_1MHz = (Vr/pi)^2 * ((65/100).^(-2) - 1)/(2*M)./burst_1MHz.k(end);

% Noise_low = ((40./100).^-2 - 1)./(2*M) .* (Vr ./ pi).^2 ./k(end);

[k_nasmyth, P_nasmyth]=univ_turb_spec(4*10^-4); %in rad/m or 1/m??
[~,ki] = min( abs(1/burst_1MHz.k(2)/2/pi - k_nasmyth) ); %NOW REMOVES .64cm of data...
ceil_1MHz = P_nasmyth(ki);
box_1MHz =  [burst_1MHz.k(2)/2/pi, burst_1MHz.k(end)/2/pi ; floor_1MHz*2*pi, ceil_1MHz]; %NOW REMOVES .64cm of data...

b2 = plot( box_1MHz(1,:), box_1MHz(2,[1,1]),'--','linewidth',2,'color',clrs(5,:) );
plot(box_1MHz(1,:),box_1MHz(2,[2,2]),'--','linewidth',2,'color',clrs(5,:))
plot(box_1MHz(1,[1,1]),box_1MHz(2,[1,2]),'--','linewidth',2,'color',clrs(5,:))
plot(box_1MHz(1,[2,2]),box_1MHz(2,[1,2]),'--','linewidth',2,'color',clrs(5,:))
plot(box_1MHz(1,:),floor2_1MHz*2*pi*[1,1],'--','linewidth',2,'color',clrs(5,:));
patch1mhz=patch([box_1MHz(1,:),fliplr(box_1MHz(1,:))],...
    [box_1MHz(2,1),box_1MHz(2,1),floor2_1MHz*2*pi,floor2_1MHz*2*pi],...
    clrs(5,:),'FaceAlpha',0.2);

grid on
title(['\epsilon = 10^{-9}, 10^{-8}, 10^{-7}, 10^{-6} [m^2 s^{-3}]'])
xlabel('Wavenumber [cpm]')
ylabel('PSD [m^2 s^{-2} cpm^{-1}]')
%legend([p1,p2,b1,b2,patch2mhz,patch1mhz],{'Pope','Nasmyth','2 MHz Range','1 MHz Range',...
%    '2 MHz Noise, $\hat{R}$ = 65-95','1 MHz Noise, $\hat{R}$ = 65-95'},...
%    'location','NorthEastOutside','interpreter','latex')
legend([p1,b1,b2,patch2mhz,patch1mhz],{'Pope','2 MHz Range','1 MHz Range',...
   '2 MHz Noise, $\hat{R}$ = 0.60-0.95','1 MHz Noise, $\hat{R}$ = 0.60-0.95'},...
   'location','NorthEastOutside','interpreter','latex') %REMOVED NASMYTH SPECTRUM IN RESPONSE TO REVIEWER COMMENTS, 4/27/21
set(gca,'ylim',10.^[-8,-2.5],'xtick',10.^[-1,0,1,2],'xlim',[1e-1,1e2])

%--------------------------------
%Epsilon Labels
% annotation('textbox',[0.2,0.5,0.1,0.1],'String','\epsilon = 1e-9','EdgeColor','None')
% ht1 = text(0.15,2e-6,'\epsilon = 10^{-9}','Rotation',-30,'Fontsize',14);
% ht2 = text(0.15,1e-5,'\epsilon = 10^{-8}','Rotation',-30,'Fontsize',14);
% ht3 = text(0.15,4.5e-5,'\epsilon = 10^{-7}','Rotation',-30,'Fontsize',14);
% ht4 = text(0.15,2.2e-4,'\epsilon = 10^{-6}','Rotation',-30,'Fontsize',14);

ht1 = text(0.11,9e-6,'\epsilon = 10^{-9} m^2/s^3','Rotation',-35,'Fontsize',14);
ht2 = text(0.11,4e-5,'\epsilon = 10^{-8} m^2/s^3','Rotation',-35,'Fontsize',14);
ht3 = text(0.11,2e-4,'\epsilon = 10^{-7} m^2/s^3','Rotation',-35,'Fontsize',14);
ht4 = text(0.11,9e-4,'\epsilon = 10^{-6} m^2/s^3','Rotation',-35,'Fontsize',14);

%Annotations 1MHz Noise Floor
% annotation('TextArrow',[0.241 0.3125],[0.3 0.35],'String','C=65 Noise',...
%     'TextBackGroundColor','w','Fontsize',13)
% annotation('Arrow',[0.241,0.28],[0.31,0.38])
% 
% annotation('Arrow',[0.31,0.28],[0.166,0.27])
% annotation('TextArrow',[0.31 0.33],[0.166 0.24],'String','C=99 Noise',...
%     'TextBackGroundColor','w','Fontsize',13)

%Annotation: Wrapping Limit
annotation('TextBox',[0.3,0.78,0.1,0.1],'String','Wrapping Limit','BackGroundColor','None',...
    'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','None','Fontsize',14)

%Annotation: Nyquist, 1/Range
text(17,1e-5,'Nyquist','FontSize',14,'Rotation',-90)
annotation('TextArrow',[0.35,0.26],[0.7,0.7],'String','[1/Range]','TextBackGroundColor','w',...
    'Fontsize',14,'VerticalAlignment','Bottom')
annotation('Arrow',[0.35,0.22],[0.73,0.73])

SPURS_figsave(gcf)

%% Unwrapping plot (current fig. 5)
clear
addpath('../processing/')

%load example burst
load('../../Data/SPURS-2/interim/NortekFiles/8116/811604_burst0510.mat')

%get histogram data
ex_index = 308;
V1_fdata= v1(ex_index,3:end);
Vmax = Vr;
SetL = length(V1_fdata(1,:));
profile = 1; %hack to make code run (most copied from histogram_unwrap_function5.m)
twopi=2*Vmax;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot and define repeating set of points%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:5
    Vrep(n,:) = V1_fdata(profile,:) + n*twopi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set proper bins and create histrograms%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins1 = twopi./50:twopi./50:twopi;
bins = [];
for nn = 1:n
    bins = [bins bins1+(nn-1).*twopi];
end
histnums = histcounts(Vrep(:),bins);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solve for largest zero window in histogram data%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Alternative to AA loops in histogram_unwrap_function3:
%
%Sometimes the biggest block of zeros is associated
% with first histogram cluster missing some values because of wrapping
%This kludge deals with that:
histnums(1:2:49)=0;histnums(2:2:50)=NaN;

ff=find(histnums~=0);
ffi=find(diff(ff)==max(diff(ff)));
win1=ff(ffi(1))+1;
win2=win1+50;

%Identify bins that were originally NaNs and ones that weren't
nan_ind=find(isnan(V1_fdata(profile,:)));
non_nan_ind=find(~isnan(V1_fdata(profile,:)));

%Place the unwrapped data in the appropriate bins
%SFZ 4/21/2020 - this section was hoarking because win2 could be longer
%than bins. Not sure what is causing this behavior...
if win2<length(bins)
    TT=find(Vrep>=bins(win1) &Vrep<bins(win2));
else %downshift 50 bins...
    win2=win2-50;
    win1=win1-50;
    TT=find(Vrep>=bins(win1) &Vrep<bins(win2));
end
%If the # of pts in TT doesn't equal # non-NaNs,
%  we have a problem.  This can be caused by amibiguous
%  identification of biggest chunk of zeros in hist
%  (e.g., two chunks same size within one wrapping range)
if length(TT)==length(non_nan_ind)
  unwrap_data(profile,non_nan_ind) = Vrep(TT)';
  unwrap_data(profile,nan_ind) = NaN;
else
  unwrap_data(profile,:) = NaN;
end
% ----------------------------



figure(5),clf
set(gcf,'Position',[114   175   917   605])
% ---------------------------- Top Row Plots
ax1=subplot(3,3,1:2);
plot(Vrep(3,:),'linewidth',1.5)
hold on
plot(Vrep([1,2,4,5],:)','-','linewidth',1.5)
grid on
plot(unwrap_data(profile,:),'ko-','linewidth',1.5,'markersize',8)
plot([0,70],[1,1]*(bins(win1)+Vr*2/3),'--k')
plot([0,70],[1,1]*(bins(win2)+Vr*2/3),'--k')
xlabel('Range Cell')
ylabel('Wrapped Velocity [m s^{-1}]')
title(['Burst0510     <T>:' datestr(mean(time)) ])

ax2=subplot(3,3,3);
histnums2 = hist(Vrep(:),bins);
plot(histnums2,bins,'linewidth',2); 
set(gca,'xlim',[0,10],'ylim',ax1.YLim)
xlabel('Counts')
ylabel('Wrapped Velocity [m s^{-1}]')
title('Histogram')

% ---------------------------- Bottom Plot

subplot(3,3,4:6)
pcolor(vUnwrap');
shading flat
hold on
plot(ex_index,63,'rv','linewidth',2,'markersize',8)
plot(ex_index,1,'r^','linewidth',2,'markersize',8)

xlabel('Burst Index')
ylabel('Range Cell')

cb = colorbar;
ylabel(cb,'Mean-removed Velocity [m s^{-1}]')
caxis([-1,1]*0.03)
title('Unwrapped')

subplot(3,3,7:9)
pcolor(v1(:,3:end)'-mean(v1(:,3:end)'));
shading flat
hold on
plot(ex_index,63,'rv','linewidth',2,'markersize',8)
plot(ex_index,1,'r^','linewidth',2,'markersize',8)

xlabel('Burst Index')
ylabel('Range Cell')
title(['Wrapped'])
cb = colorbar;
ylabel(cb,'Mean-removed Velocity [m s^{-1}]')
caxis([-1,1]*0.03)

% % ---------------------------- Inset Plot
% ax_inset = axes('Position',[.2 .13 .1 .38]);
% plot(Vrep(3,:),1:63,'linewidth',1.5)
% hold on
% plot(Vrep([1,2,4,5],:)',1:63,'-','linewidth',1.5)
% grid on
% plot(unwrap_data(profile,:),1:63,'ko-','linewidth',1.5,'markersize',8)
% set(gca,'xlim',[0.29,0.335],'ylim',[0,64],'XaxisLocation','top')

x_left = 0.13;
y_top = 0.9;

annotation('Textbox',[x_left,y_top,0.015,0.025],'String','a',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
annotation('Textbox',[0.695,y_top,0.015,0.025],'String','b',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
annotation('Textbox',[x_left,0.6,0.015,0.025],'String','c',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
annotation('Textbox',[x_left,0.3,0.015,0.025],'String','d',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')

SPURS_figsave(gcf)

%% Pot. Flow Plot (current fig. 6)

clear
close all

clrs = colororder;

%ts file for "instrument" data
%=load('../../Data/SPURS-2/interim/NortekFiles/8089/8089_dissipation_timeseries.mat'); %old burst
ts =load('../../Data/SPURS-2/interim/NortekFiles/8116/8116_dissipation_timeseries.mat');
blanking = ts.instrument.range(1)-ts.instrument.binSize;
BinSz = ts.instrument.binSize;

for bi = 1000
    %OLD - wanted same burst as unwrapping example!
    %fn = ['808903_burst' num2str(bi,'%04.f') '.mat'];
    %burst = load(['../../Data/SPURS-2/interim/NortekFiles/8089/' fn]);
    
    burst = load('../../Data/SPURS-2/interim/NortekFiles/8116/811604_burst0510.mat');


%Model Parrams OLD, for 808903_burst1000.mat. Changed to match burst of
%Fig. 5

% R = 0.3;
% U_inf = 0.1;
% r = linspace(R,2.5,20);
% inds = 215:230;

R = 0.3;
U_inf = 0.1;
r = linspace(R,2.5,20);
inds = 245:265;

figure(6),clf
hold on
% ln1=plot(ts.instrument.range+R, burst.vUnwrap(215:230,:)');
% model = U_inf.*(1- (R)^2./(ts.instrument.range+R).^2 );
% model_extended = U_inf.*(1- (R)^2./(r).^2 );
% ln2=plot(ts.instrument.range+R, mean(model)-model,'k','linewidth',2);
% plot(r, mean(model)-model_extended,'k--','linewidth',2);
% ln3 = plot([1,1]*(0.65+R),[-0.05,0.1],'-r','linewidth',2);
% ln4 = plot([1,1]*R,[-0.05,0.1],'-','linewidth',2);
% set(gca,'ylim',[-0.05,0.1])

model = U_inf.*(1- (R)^2./(ts.instrument.range+R).^2 );
model_extended = U_inf.*(1- (R)^2./(r).^2 );

ln1=plot(ts.instrument.range+R, burst.vUnwrap(inds,:)' - mean(model));
ln2=plot(ts.instrument.range+R, -model,'k','linewidth',2);
plot(r, -model_extended,'k--','linewidth',2);

ln3 = plot([1,1]*(0.65+R),[-0.15,0],'-r','linewidth',2);
ln4 = plot([1,1]*R,[-0.15,0],'-','linewidth',2,'color',clrs(6,:));
set(gca,'ylim',[-U_inf*1.25,0])


grid on
legend([ln1(1),ln2,ln4,ln3],{'Meas.','Pot. Flow Model','AQD Location','V/V_\infty = 0.9'})
xlabel('Range from mooring line [m]')
ylabel('Velocity [m/s]')
title('Burst0510, SPURS-2, 41.5m Depth, index: 245:265')
end

SPURS_figsave(gcf)

%% Spectral Fits (current fig. 7)
clear
close all
addpath('../processing/')
%burst = load('../../Data/SPURS-2/interim/NortekFiles/8116/811604_burst0510.mat');
%burst = load('../../Data/SPURS-1/interim/NortekFiles/5143/c_000001_burst1020.mat');
%load('/Users/zippelsf/Documents/Projects/SPURS/Data/SPURS-1/interim/NortekFiles/5143/5143_dissipation_timeseries.mat','instrument');

%CHANGED TO SAME BURST AS OTHER EXAMPLES
burst = load('../../Data/SPURS-1/interim/NortekFiles/6774/677404_burst1865.mat');
load('/Users/zippelsf/Documents/Projects/SPURS/Data/SPURS-1/interim/NortekFiles/6774/6774_dissipation_timeseries.mat','instrument');

%burst = load('../../Data/SPURS-2/interim/NortekFiles/8116/811604_burst0510.mat');
%load('../../Data/SPURS-2/interim/NortekFiles/8116/8116_dissipation_timeseries.mat','instrument');

%Burst Fit
L1 = instrument.Lxmit;
L2 = instrument.Lrecv;

%Burst_fit = ( (18/55)*(8/9/0.4)^(2/3).*burst.epsilon^(2/3).*burst.k.^(-5/3) + burst.N(1) ) ...
%            .*sinc(L1.*burst.k./2./pi).^2 .* sinc(L2.*burst.k./2./pi).^2 + burst.N(2);
Burst_fit = ( (18/55)*(8/9/0.4)^(2/3).*burst.epsilon^(2/3).*burst.k.^(-5/3) + burst.N(1) ) ...
            .*sinc(L1.*burst.k./2./pi).^2 .* sinc(L2.*burst.k./2./pi).^2;

figure(7),clf
set(gcf,'Position',[150   368   955   384])
subplot(1,2,1)
loglog(burst.k(2:end-1)/2/pi, burst.Pxx(2:end-1)*2*pi,'linewidth',2)
hold on
plot(burst.k(2:end)/2/pi, Burst_fit(2:end)*2*pi,'-k','linewidth',2)

legend({'Measured Spectrum','Fit Spectrum'})
grid on
xlabel('Wavenumber [cpm]')
ylabel('PSD [m^2 s^{-2} cpm^{-1}]')
title('Example Burst1865, SPURS-1, 21.5m')
set(gca,'ylim',10.^[-6.8,-4.3])

subplot(1,2,2)
loglog(burst.k(2:end-1)/2/pi, (burst.Pxx(2:end-1))*2*pi ./sinc(L1.*burst.k(2:end-1)./2./pi).^2 ./ sinc(L2.*burst.k(2:end-1)./2./pi).^2,'linewidth',2)
hold on
plot(burst.k(2:end)/2/pi, ( (18/55)*(8/9/0.4)^(2/3).*burst.epsilon^(2/3).*burst.k(2:end).^(-5/3) + burst.N(1) )*2*pi,'-k','linewidth',2)
plot([burst.k(2),burst.k(end) ]/2/pi, [1,1]*burst.N(1)*2*pi,'--r','linewidth',2)
plot(burst.k/2/pi, ( (18/55)*(8/9/0.4)^(2/3).*burst.epsilon^(2/3).*burst.k.^(-5/3))*2*pi,'--','linewidth',2 )

legend({'Corrected Spectrum','Fit Spectrum','White Noise Estimate','Intertial Subrange Estimate'})
grid on
xlabel('Wavenumber [cpm]')
ylabel('PSD [m^2 s^{-2} cpm^{-1}]')
%title('Example Burst0510, SPURS-2, 7m')
set(gca,'ylim',10.^[-6.8,-4.3])

x_left = 0.13;
y_top = 0.875;
annotation('Textbox',[x_left,y_top,0.02,0.05],'String','a',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[0.57,y_top,0.02,0.05],'String','b',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')

SPURS_figsave(gcf)

%% Fig. 8 are microstructure comparisons /w. tplots ... see "compare_with_glider.mlx"
clear
close all

%load SPURS-1 files
PROJ = 'SPURS-1';
fp = ['../../Data/' PROJ '/interim/NortekFiles/'];
load([fp 'SPURS1_dissipation_grid_v1c.mat']); %NOTE
load(['../../Data/' PROJ '/interim/BuoyancyFlux_b.mat']);
met = load(['../../Data/' PROJ '/processed/spurs1_met_1hr.mat'],'sst','sal','mday');

z_L = (dissipation.depth ./ L_pen')';
clrs = get( groot, 'defaultAxescolororder'); %COLORS
depth_clrs = [[1,1,1]*0.2 ;[1,1,1]*0.6];

%Load Gliderr Data
glider = load('../../Data/SPURS-1/external/StLaurent_turbulence/epC2_grid.mat');
fp = '/Users/zippelsf/Documents/Projects/SPURS/Data/SPURS-1/external/StLaurent_turbulence/tem1_grid.mat';
glider_temp = load(fp);

%Load SPURS-1 Temperature Data
load('/Users/zippelsf/Documents/Projects/SPURS/Data/SPURS-2/processed/SPURS1_TS_grid.mat','S');
%S.TIME is days since 1950-01-01T00:00:00Z
time_matlab_ref = S.TIME+datenum('1950-01-01 00:00:00');
inds = time_matlab_ref>glider_temp.timvec(1) & time_matlab_ref<glider_temp.timvec(end);
inds2 = dissipation.time>glider_temp.timvec(1) & dissipation.time<glider_temp.timvec(end);
neg_bflux_inds=dissipation.time>glider_temp.timvec(1) & dissipation.time<glider_temp.timvec(end) & L'<0;
mooring_inds = dissipation.time > min(glider.timvec) & dissipation.time < max(glider.timvec);
mooring_temp_inds = S.TIME+datenum('1950-01-01 00:00:00') > min(glider.timvec) & S.TIME+datenum('1950-01-01 00:00:00') < max(glider.timvec);

%----------------MAKE THIS FASTER BY ONLY INCLUDING RELEVANT INDEXING
sst_interp = interp1(met.mday, met.sst, S.TIME(mooring_temp_inds)+datenum('1950-01-01 00:00:00'));
sss_interp = interp1(met.mday, met.sal, S.TIME(mooring_temp_inds)+datenum('1950-01-01 00:00:00'));

T_interp = interp1(cat(1,0,S.DEPTH), cat(1,sst_interp' ,S.TEMP(:,mooring_temp_inds)), 0:1:250)';
S_interp = interp1(cat(1,0,S.DEPTH), cat(1,sss_interp' ,S.PSAL(:,mooring_temp_inds)), 0:1:250)';


c_range = [26.9,27.8];

figure(8),clf
set(gcf,'Position',[-1527         938        1324         629])

ax1=subplot(4,1,1);
[h,c]=contourf(glider_temp.timvec, glider_temp.pressvec, glider_temp.propmat+6.6368,linspace(24,29,100));
set(c,'LineColor','none');
axis ij
cb=colorbar;
ylabel(cb,'Temperature')
ylabel('Depth')
caxis(c_range)
hold on
plot([datenum('09/23/2012'),datenum('10/5/2012')], [1,1]*12.5,'--','linewidth',2,'color',depth_clrs(1,:))
plot([datenum('09/23/2012'),datenum('10/5/2012')], [1,1]*21.5,'--','linewidth',2,'color',depth_clrs(2,:))
%plot(dissipation.time(inds2), ml_depth_04(inds2),'-k','linewidth',2.5)
%plot(dissipation.time(inds2), ml_depth_04(inds2),'-w','linewidth',1)
%plot(dissipation.time(neg_bflux_inds), dissipation.time(neg_bflux_inds)*0+1,'*r','linewidth',1)
set(gca,'xlim',[datenum('9/23/2012'),datenum('10/04/2012 12:00:00')],'ylim',[0,30],'xtick',xticks)
title('Glider')
datetick('x','keeplimits','keepticks')

ax2=subplot(4,1,2);
[h,c]=contourf(time_matlab_ref(inds), (0:1:100), T_interp(:,1:101)',linspace(24,29,100));
set(c,'LineColor','none');
axis ij
cb=colorbar;
ylabel(cb,'Temperature')
ylabel('Depth')
caxis(c_range)

hold on
plot([datenum('09/23/2012'),datenum('10/5/2012')], [1,1]*12.5,'--','linewidth',2,'Color',depth_clrs(1,:))
plot([datenum('09/23/2012'),datenum('10/5/2012')], [1,1]*21.5,'--','linewidth',2,'color',depth_clrs(2,:))
%plot(dissipation.time(inds2), ml_depth_04(inds2),'-k','linewidth',2.5)
%plot(dissipation.time(inds2), ml_depth_04(inds2),'-w','linewidth',1)
%plot(dissipation.time(neg_bflux_inds), dissipation.time(neg_bflux_inds)*0+1,'*r','linewidth',1)
set(gca,'xlim',[datenum('9/23/2012'),datenum('10/04/2012 12:00:00')],'ylim',[0,30],'xtick',xticks)
title('Mooring')
datetick('x','keeplimits','keepticks')

% % %----------------------- TDIFF SUBPLOT
% % subplot(5,1,3)
% % for di = 1:2%:length(dissipation.depth)
% %     hold on
% %     glider_inds = glider_temp.pressvec == (dissipation.depth(di));
% %     %interp to same grid
% %     dt_mooring_interp = interp1(time_matlab_ref(inds), T_interp(:,4)'-T_interp(:,round(dissipation.depth(di)))', glider_temp.timvec);
% %     dt_glider = nanmean(glider_temp.propmat(7,:),1)-glider_temp.propmat(glider_inds,:);
% %     
% %     %Only look at "stable" profiles 
% %     %dt_mooring_interp(dt_mooring_interp<0) = NaN;
% %     %dt_glider(dt_glider<0) = NaN;
% %     
% %     plot(glider_temp.timvec, movmean((dt_glider-dt_mooring_interp'),11,'omitnan'),'linewidth',2,'color',depth_clrs(di,:))
% % 
% %     set(gca,'xlim',[datenum('9/23/2012'),datenum('10/04/2012 12:00:00')])
% %     datetick('x','keeplimits')
% %     grid on
% %     ylabel('\Delta T_{gldr} - \Delta T_{mring}')
% % end
% % legend({'12.5m','21.5m'},'location','best')
% % set(gca,'xlim',[datenum('9/23/2012'),datenum('10/04/2012 12:00:00')])
% % grid on
% % cb1 = colorbar;
% % set(cb1,'Visible','off')
% % title('Difference in Temperature Gradient [from 3m ref]')

%----------------------- Dissipation SUBPLOTs
t_offset = -15/60/24; %guessing here
for di = 1:2%:length(dissipation.depth)
    subplot(4,1,di+2)
    hold on
    glider_inds = glider.pressvec == (dissipation.depth(di));
    mooring_inds = dissipation.time > min(glider.timvec) & dissipation.time < max(glider.timvec) & dissipation.Spectra(di).DOF'>9;
    
    zL_inds = mooring_inds & z_L(di,:)<-10 & (dissipation.depth(di)<0.75*ml_depth_04');
    if  ~isempty(zL_inds)
        ln_jb=plot(dissipation.time(zL_inds)+t_offset,Jb(zL_inds),'s','MarkerSize',5,'LineWidth',1.5,'color',[1,1,1]*0.5,'MarkerFaceColor',[1,1,1]*0.5);
    end
    
    gldr=semilogy(glider.timvec, movmean(glider.propmat(glider_inds,:),8,'omitnan'),'-x','color',clrs(1,:),'linewidth',2);
    mring=semilogy(dissipation.time(mooring_inds)+t_offset, (dissipation.epsilon(di,mooring_inds)),'-o','linewidth',2,'color',clrs(2,:));
    
    %Grey box below 1e-8
    Xvertex = [datenum('9/23/2012'),datenum('9/23/2012'),datenum('10/04/2012 12:00:00'),datenum('10/04/2012 12:00:00')];
    Yvertex = 10.^[-8,-10,-10,-8];
    fill(Xvertex,Yvertex,[1,1,1]*0.6,'facealpha',0.5,'edgecolor','none');
    
    title(['Depth: ' num2str(dissipation.depth(di),3)])
    set(gca,'yscale','log','xlim',[datenum('9/23/2012'),datenum('10/04/2012 12:00:00')])
    datetick('x','keeplimits')
    leg = legend([gldr,mring,ln_jb],{'Glider','Mooring','Jb_0'});
    leg.Position = [0.87 0.1739 0.0616 0.0771];
    grid on
    ylabel('\epsilon [m^2/s^3]')
    cb1 = colorbar;
    set(cb1,'Visible','off')
end

cmap = parula(22);
%cmap = cptcmap('thermal');
colormap(cmap)

% ANNOTATIONS
%figure labels
x_left = 0.135;
annotation('Textbox',[x_left,0.89,0.014,0.03],'String','a',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[x_left,0.67,0.014,0.03],'String','b',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[x_left,0.45,0.014,0.03],'String','c',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[x_left,0.24,0.014,0.03],'String','d',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[x_left,0.24,0.014,0.03],'String','e',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')

%boxes
box_a1 = annotation('textbox',[0.235, 0.8, 0.03, 0.08],'string','A.1','FitBoxToText','off');
box_b1 = annotation('textbox',[0.235, 0.58, 0.03, 0.08],'string','B.1','FitBoxToText','off');
box_d1 = annotation('textbox',[0.235, 0.115, 0.03, 0.15],'string','D.1','FitBoxToText','off');

%9/27 box
box_a2 = annotation('textbox',[0.415, 0.8, 0.03, 0.08],'string','A.2','FitBoxToText','off');
box_b2 = annotation('textbox',[0.415, 0.58, 0.03, 0.08],'string','B.2','FitBoxToText','off');
box_c2 = annotation('textbox',[0.415, 0.34, 0.03, 0.13],'string','C.2','FitBoxToText','off');

%9/29 box
box_a3 = annotation('textbox',[0.5, 0.8, 0.04, 0.1],'string','A.3','FitBoxToText','off');
box_b3 = annotation('textbox',[0.5, 0.58, 0.04, 0.1],'string','B.3','FitBoxToText','off');
box_c3 = annotation('textbox',[0.5, 0.34, 0.04, 0.13],'string','C.3','FitBoxToText','off');
box_d3 = annotation('textbox',[0.5, 0.115, 0.04, 0.12],'string','D.3','FitBoxToText','off');

%9/30 box
box_a4 = annotation('textbox',[0.57, 0.8, 0.04, 0.14],'string','A.4','FitBoxToText','off');
box_b4 = annotation('textbox',[0.57, 0.58, 0.04, 0.14],'string','B.4','FitBoxToText','off');
box_c4 = annotation('textbox',[0.57, 0.34, 0.04, 0.14],'string','C.4','FitBoxToText','off');
box_d4 = annotation('textbox',[0.57, 0.115, 0.04, 0.14],'string','D.4','FitBoxToText','off');

%Box at end?
box_a5 = annotation('textbox',[0.73, 0.78, 0.13, 0.16],'string','A.5','FitBoxToText','off');
box_b5 = annotation('textbox',[0.73, 0.56, 0.13, 0.16],'string','B.5','FitBoxToText','off');
box_c5 = annotation('textbox',[0.73, 0.34, 0.13, 0.14],'string','C.5','FitBoxToText','off');
box_d5 = annotation('textbox',[0.73, 0.115, 0.13, 0.14],'string','D.5','FitBoxToText','off');

SPURS_figsave(gcf)

%% Glider mooring Compare [removed in place of scatterplot]
%load SPURS-1 files
PROJ = 'SPURS-1';
fp = ['../../Data/' PROJ '/interim/NortekFiles/'];
load([fp 'SPURS1_dissipation_grid_v1c.mat']); %NOTE v1b used in JTECH
load(['../../Data/' PROJ '/interim/BuoyancyFlux_b.mat']);
z_L = (dissipation.depth ./ L_pen')';

glider = load('../../Data/SPURS-1/external/StLaurent_turbulence/epC2_grid.mat');
fp = '/Users/zippelsf/Documents/Projects/SPURS/Data/SPURS-1/external/StLaurent_turbulence/tem1_grid.mat';
glider_temp = load(fp);

xticks = [floor(glider_temp.timvec(1)):1:ceil(glider_temp.timvec(end))];
 clrs = colororder;
mooring_inds = dissipation.time > min(glider.timvec) & dissipation.time < max(glider.timvec);

%---------------------------------------- GET GLIDERR DISTANCE TO MOORRIING
load('../../Data/SPURS-1/external/StLaurent_turbulence/helo_data.mat','g')
load('../../Data/SPURS-1/external/StLaurent_turbulence/SPURS1_metob.mat')

t1 = datenum('9/22/2012');
t2 = datenum('10/5/2012');
gi = find(g.m_present_time > t1 & g.m_present_time < t2);
bi = find(metob.mday > t1 & metob.mday < t2);

g_lon = dm2degrees( [fix(g.m_lon(gi)/100), abs(-(g.m_lon(gi)/100) + fix(g.m_lon(gi)/100))*100] );
g_lat = dm2degrees( [floor(g.m_lat(gi)/100), abs((g.m_lat(gi)/100) - fix(g.m_lat(gi)/100))*100] );
% 
% clf
% plot(g_lon, g_lat,'o')
% hold on
% plot(metob.GPSLON(bi), metob.GPSLAT(bi),'x')
% legend({'Glider','Mooring'},'location','best')
% grid on

buoy_ll_interp = interp1(metob.mday(bi),[metob.GPSLAT(bi); metob.GPSLON(bi)]',g.m_present_time(gi));
earthRadiusInMeters = 6371000;
gldr_mooring_dist_m =  distance(buoy_ll_interp(:,1), buoy_ll_interp(:,2), g_lat, g_lon, earthRadiusInMeters);

gldr_dist_interp = interp1(g.m_present_time(gi),gldr_mooring_dist_m, glider.timvec);
mooring_dist_interp = interp1(g.m_present_time(gi),gldr_mooring_dist_m, dissipation.time);

clear g metobs
%----------------------------------------
glider_timemat = (glider.timvec*ones(1,length(glider.pressvec)))';
glider_epsilon_interp = nan(500,4);
mooring_epsilon_small_mat = nan(500,4);
mooring_small_tvec = nan(500,4);
max_dist = 6000;

figure(8)
clf
set(gcf,'Position',[ -1674,892,1248,753])

t_offset = -(15/60)/24; %dissipation has been shifted roughly +15min to match flux files
for di = 1:4%:length(dissipation.depth)
    ax(di) = subplot(4,1,di);
    hold on
    glider_inds = abs(glider.pressvec - dissipation.depth(di))<1 & gldr_dist_interp'<max_dist;
    %glider_inds = abs(glider.pressvec - dissipation.depth(di))<0.5; %glider vec index (instead of mat index)
    %low_quality = (dissipation.epsilon(di,:)').^(2/3) < ...
    %((18/55).*(8/9/0.4)^(2/3).*dissipation.Spectra(di).k(:,3).^(-5/3)).^(-1).*dissipation.Spectra(di).N(:,1)./sqrt(19).*dissipation.Spectra(di).DOF(:) ./ chi2inv(0.05,dissipation.Spectra(di).DOF(:));
    %mooring_inds = dissipation.time > min(glider.timvec) & dissipation.time < max(glider.timvec) &...
    %   (~low_quality)';
    mooring_inds = dissipation.time > min(glider.timvec) & dissipation.time < max(glider.timvec) &...
        (dissipation.Spectra(di).DOF>9)';
    
    %ln_jb = plot(dissipation.time(mooring_inds)+t_offset, Jb(mooring_inds), '-','LineWidth',3,'color',[1,1,1]*0.3); %FULL JB LINE 
    ln_gldr = plot(glider_timemat(glider_inds), movmean(glider.propmat(glider_inds),8,'omitnan'),'x','color',clrs(1,:),'linewidth',2);
    %ln_gldr = plot(glider.timvec, movmean(glider.propmat(glider_inds,:),11,'omitnan'),'-x','color',clrs(1,:),'linewidth',2);%OLD PLOTTING 
    
    ln_mrng = semilogy(dissipation.time(mooring_inds)+t_offset, dissipation.epsilon(di,mooring_inds),'-o','markerfacecolor',clrs(2,:),'color',clrs(2,:));   
    
    zL_inds = mooring_inds & z_L(di,:)<-10 & (dissipation.depth(di)<0.9*ml_depth_04');
    if  ~isempty(zL_inds)
        ln_zl = plot(dissipation.time(zL_inds), Jb(zL_inds),'o','MarkerSize',8,'LineWidth',2,'color',clrs(3,:));
    end
    
    title(['Depth: ' num2str(dissipation.depth(di),3)])
    set(gca,'yscale','log')
    datetick('x')
    %legend([ln_gldr,ln_mrng,ln_jb,ln_zl],{'Glider','Mooring','Jb_0','Jb_0 (z/L < -5)'},'location','EastOutside','Fontsize',14)
    legend([ln_gldr(1),ln_mrng,ln_zl],{'Glider','Mooring','Jb_0'},'location','EastOutside','Fontsize',14)
    grid on
    ylabel('\epsilon [m^2/s^3]')
    
    %-------------------------- TRY PLOTTING DIRECT scatter
    %avg all depths at same time, then all times?
    gldr_depths_i = find(abs(glider.pressvec - dissipation.depth(di))<1);
    gldr_times_i = find(gldr_dist_interp'<max_dist);
    glider_epsilon_interp(1:sum(mooring_inds),di) = interp1(glider.timvec(gldr_times_i), movmean(nanmean(glider.propmat(gldr_depths_i,gldr_times_i),1),11,'omitnan'), dissipation.time(mooring_inds)+t_offset,'linear')';
    mooring_epsilon_small_mat(1:sum(mooring_inds),di) = dissipation.epsilon(di,mooring_inds)';
    mooring_small_tvec(1:sum(mooring_inds),di) = dissipation.time(mooring_inds)+t_offset;
end

linkaxes(ax,'x')
%Annotations
x_left = 0.13;
annotation('Textbox',[x_left,0.9,0.015,0.02],'String','a',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[x_left,0.68,0.015,0.02],'String','b',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[x_left,0.46,0.015,0.02],'String','c',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[x_left,0.25,0.015,0.02],'String','d',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
%SPURS_figsave(gcf)

%% Mooring-glider Scatterplot (Run Previous section first w/o clearing)

%  figure(100),clf

concat_glider=[];
concat_mooring = [];
for di = 1:3
% subplot(1,4,di)
% hold on
%Remove suspect indexes
if di==1
bi = (mooring_small_tvec(:,di) > datenum('10/2/2012 06:00')) |...
     ( mooring_small_tvec(:,di) > datenum('9/28/2012 19:00') & mooring_small_tvec(:,di) < datenum('9/29/2012 06:00'))|...
     ( mooring_small_tvec(:,di) > datenum('9/29/2012 21:00') & mooring_small_tvec(:,di) < datenum('9/30/2012 10:00'));
elseif di==2
bi = (mooring_small_tvec(:,di) > datenum('10/2/2012 06:00')) |...
     ( mooring_small_tvec(:,di) > datenum('9/28/2012 19:00') & mooring_small_tvec(:,di) < datenum('9/29/2012 06:00')) | ...
     ( mooring_small_tvec(:,di) > datenum('9/29/2012 21:00') & mooring_small_tvec(:,di) < datenum('9/30/2012 10:00')) |...
     ( mooring_small_tvec(:,di) > datenum('9/27/2012 21:00') & mooring_small_tvec(:,di) < datenum('9/28/2012 10:00')) |...
     ( mooring_small_tvec(:,di) > datenum('9/24/2012 17:00') & mooring_small_tvec(:,di) < datenum('9/25/2012 2:00'));
elseif di>2
bi = (mooring_small_tvec(:,di) > datenum('10/2/2012 06:00')) |...
     (mooring_small_tvec(:,di) < datenum('9/27/2012 12:00'));
end
% loglog(glider_epsilon_interp(bi,di),mooring_epsilon_small_mat(bi,di),'x','color',[1,1,1]*0.7)
% loglog(glider_epsilon_interp(~bi,di),mooring_epsilon_small_mat(~bi,di),'.','color',[1,1,1]*0.3)
% [Xbin,Ybin,Ylo,Yhi,CI] = geo_bin_data(glider_epsilon_interp(~bi,di), mooring_epsilon_small_mat(~bi,di), logspace(-9,-6,18), 10);
% 
% plot(Xbin,Ybin,'d','markersize',10,'MarkerFaceColor',clrs(di,:),'color',clrs(di,:))
% plot([Xbin;Xbin],[Ybin./CI;Ybin.*CI], '-','LineWidth',1,'color',clrs(di,:))
% 
% plot(10.^[-10,-6],10.^[-10,-6],'--k','linewidth',2)
% xlabel('Glider \epsilon [m^2 s^{-3}]')
% ylabel('Mooring \epsilon [m^2 s^{-3}]')
% set(gca,'yscale','log','xscale','log')
% title(['Depth: ' num2str(dissipation.depth(di),3)])

if di < 3
concat_glider = cat(1,concat_glider, glider_epsilon_interp(~bi,di));
concat_mooring = cat(1,concat_mooring, mooring_epsilon_small_mat(~bi,di));
end

end

%BIN ALL 12m & 21m DATA TOGETHER ---------
figure(9),clf
hold on

plot(10.^[-10,-6],2.*10.^[-10,-6],'--','linewidth',2,'color',[1,1,1]*0.7)
plot(10.^[-10,-6],0.5.*10.^[-10,-6],'--','linewidth',2,'color',[1,1,1]*0.7)

% loglog(concat_glider, concat_mooring,'.','color',[1,1,1]*0.3)
%[Xbin,Ybin,Ylo,Yhi,CI] = geo_bin_data(concat_glider, concat_mooring, logspace(-9,-6,20), 10);
[Xbin,Ybin,Ylo,Yhi,CI] = geo_bin_data(concat_glider, concat_mooring, logspace(-10,-6,40), 10);

plot(Xbin,Ybin,'d','markersize',10,'MarkerFaceColor','k','color','k')
%plot([Xbin;Xbin],[Yhi;Ylo], '-','LineWidth',1,'color',clrs(di,:))
plot([Xbin;Xbin],[Ybin./CI;Ybin.*CI], '-','LineWidth',2,'color','k')

grid on
plot(10.^[-10,-6],10.^[-10,-6],'--k','linewidth',2)

xlabel('Glider \epsilon [m^2 s^{-3}]')
ylabel('Mooring \epsilon [m^2 s^{-3}]')
set(gca,'yscale','log','xscale','log','xlim',[1e-10,1e-6],'ylim',[1e-10,1e-6])
title(['Combined Depths 12.5, 21.5, (suspect data excluded)'])

SPURS_figsave(gcf)

%% Figures 10 and 11 - 
clear
close all

PROJ = 'SPURS-1';

fp = ['../../Data/' PROJ '/interim/NortekFiles/'];
load([fp 'SPURS1_dissipation_grid_v1c.mat']);
%old=load([fp 'SPURS1_dissipation_grid_v1b.mat']);
load(['../../Data/' PROJ '/interim/BuoyancyFlux_c.mat']);
flux = load(['../../Data/' PROJ '/processed/spurs1_flux_1hr.mat'],'Qr');
met = load(['../../Data/' PROJ '/processed/spurs1_met_1hr.mat']);


z_L = (dissipation.depth ./ L_pen')';
%z_L = (dissipation.depth ./ L)';
clrs = colororder;

%load Glider
gldr = load('../../Data/SPURS-1/external/StLaurent_turbulence/all.mat');

mooring_inds = dissipation.time > min(min(gldr.DAT)) & dissipation.time < max(max(gldr.DAT));
figure(10)
clf
set(gcf,'Position',[16,1153,953,655]);
%mean_gldr_interp = nan(7,500)
ii = 0;
for di = [1,2,3,4,5,7]%:length(dissipation.depth)
    ii = ii+1;
    subplot(6,1,ii)
    hold on
    
    glider_inds = find(abs(gldr.PRS - dissipation.depth(di)) < 1.5);
    zl_inds = find(mooring_inds & z_L(di,:)<-5 & dissipation.depth(di)< ml_depth_04');
    
    semilogy(gldr.DAT(glider_inds), movmean(gldr.EPS_qc(glider_inds),8,'omitnan'),'-x','linewidth',1.5,'color',clrs(1,:))
    %semilogy(dissipation.time(mooring_inds), dissipation.epsilon(di,mooring_inds),'-o','linewidth',1.5,'color',clrs(2,:))
    semilogy(dissipation.time(mooring_inds), dissipation.epsilon(di,mooring_inds),'-o','linewidth',1.5,'color',clrs(2,:))
    
    title(['Depth: ' num2str(dissipation.depth(di),3)])
    set(gca,'yscale','log')
    datetick('x')
    ylabel('\epsilon [m^2/s^3]')
    %legend({'Glider','Mooring','Jb_0'},'location','best')
    legend({'Glider','Mooring'},'location','East Outside')
    grid on
end


%Interpolation for direct comparison
%Make sure not to compare interpolated values over gaps
mooring_inds = dissipation.time > min(min(gldr.DAT)) & dissipation.time < max(max(gldr.DAT));

for di = 1:7%:length(dissipation.depth)
    glider_inds = find(abs(gldr.PRS - dissipation.depth(di)) < 1);
    %make comparison index
    k_ind{di} = dsearchn(gldr.DAT(glider_inds), dissipation.time(mooring_inds)');
    mean_gldr_ep = movmean(gldr.EPS_qc(glider_inds),9,'omitnan');
    mean_gldr_interp{di} = mean_gldr_ep(k_ind{di});
    mean_gldr_interp{di}( abs(gldr.DAT(glider_inds(k_ind{di})) - dissipation.time(mooring_inds)') > 0.25/24) = NaN;
    
    
    %moring_interp{di} = dissipation.epsilon(di,mooring_inds)';
    moring_interp{di} = movmean(dissipation.epsilon(di,mooring_inds)',3,'omitnan');
    
end

SPURS_figsave(gcf)

%% Fig 11 - run w/o clearing last section
figure(11),clf

%loglog([mean_gldr_interp{[1:3,5,7]}],[moring_interp{[1:3,5,7]}],'x')
hold on
X = [mean_gldr_interp{[1:4,5,7]}];
Y = [moring_interp{[1:4,5,7]}];
[binX, binY,Ylo,Yhi,CI] = geo_bin_data(X(:),Y(:), logspace(-10,-6,40), 10);
plot(binX,binY, 'd','markersize',8,'LineWidth',2,'markerfacecolor','k','color','k')
plot([binX;binX],[binY./CI;binY.*CI], '-','LineWidth',1,'color','k')

plot(10.^[-10,-6],2.*10.^[-10,-6],'--','linewidth',2,'color',[1,1,1]*0.7)
plot(10.^[-10,-6],0.5.*10.^[-10,-6],'--','linewidth',2,'color',[1,1,1]*0.7)
plot([1e-10,1e-6],[1e-10,1e-6],'--k','linewidth',2)

grid on
set(gca,'yscale','log','xscale','log','xlim',[1e-10,1e-6],'ylim',[1e-10,1e-6])
title('SPURS-1, Glider Deployment 2')
xlabel('\epsilon Glider [m^2 s^{-3}]')
ylabel('\epsilon Mooring [m^2 s^{-3}]')

SPURS_figsave(gcf)

%% Jb0 comparison (current fig. 12)

clear
prj = {'SPURS-1','SPURS-2'};

figure(12),clf
set(gcf,'Position',[-1318        1105         639         364])

for pi = 1:2
    PROJ = prj{pi};
    di = 1; %depth index
if strcmp(PROJ,'SPURS-2')
    fp = ['../../Data/' PROJ '/interim/NortekFiles/'];
    load([fp 'SPURS2_dissipation_grid_v1c.mat']); %note, v1b used for JTECH subm.
    load(['../../Data/' PROJ '/interim/BuoyancyFlux_b.mat']);
    flux = load(['../../Data/' PROJ '/processed/spurs2_flux_1hr.mat'],'Qr');
    %dissipation.depth = dissipation.depth(1:4);
elseif strcmp(PROJ,'SPURS-1')
    fp = ['../../Data/' PROJ '/interim/NortekFiles/'];
    load([fp 'SPURS1_dissipation_grid_v1c.mat']); %note, v1b used for JTECH subm.
    load(['../../Data/' PROJ '/interim/BuoyancyFlux_c.mat']);
    flux = load(['../../Data/' PROJ '/processed/spurs1_flux_1hr.mat'],'Qr');
end
z_L = (dissipation.depth ./ L_pen')';
clrs = get( groot, 'defaultAxescolororder'); %COLORS

subplot(1,2,pi)
gi = find(z_L(di,:) < -5 & dissipation.depth(di)< 0.75*(ml_depth_04'));
hold on
plot(Jb(gi), dissipation.epsilon(di,gi),'x','color',[1,1,1]*0.6)
[binX, binY,Ylo,Yhi,CI] = geo_bin_data(Jb(gi), dissipation.epsilon(di,gi)', logspace(-8,-6,21), 5);
plot(binX,binY, 'd','markersize',8,'LineWidth',2,'markerfacecolor','k','color','k')
plot([binX;binX],[Ylo;Yhi], '-','LineWidth',1,'color','k')

xlabel('Jb_0 [m^2/s^3]')
ylabel('\epsilon [m^2/s^3]')
set(gca,'yscale','log','xscale','log','xlim',[3e-9,5e-7])
title([PROJ ' Depth: ' num2str(dissipation.depth(1)) ])
grid on
hold on
plot(10.^[-9,-6],10.^[-9,-6],'--k','linewidth',2)
end

x_left = 0.135;
y_top = 0.875;
annotation('Textbox',[x_left,y_top,0.02,0.05],'String','a',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')
annotation('Textbox',[0.575,y_top,0.02,0.05],'String','b',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',18,'BackgroundColor','w')

SPURS_figsave(gcf)

%% depth-time plot w/ zoom (current fig. 13)

clear

PROJ = 'SPURS-1';
figure(13),clf
set(gcf,'Position',[80,75,1187,730])

if strcmp(PROJ,'SPURS-2')
    fp = ['../../Data/' PROJ '/interim/NortekFiles/'];
    load([fp 'SPURS2_dissipation_grid_v1c.mat']);
    load(['../../Data/' PROJ '/interim/BuoyancyFlux_b.mat']);
    %flux = load(['../Data/' PROJ '/processed/spurs2_flux_1hr.mat'],'Qr');
    %dissipation.depth = dissipation.depth(1:4);
elseif strcmp(PROJ,'SPURS-1')
    fp = ['../../Data/' PROJ '/interim/NortekFiles/'];
    load([fp 'SPURS1_dissipation_grid_v1c.mat']); %NOTEE v1b was used in JTECH initial sub.
    load(['../../Data/' PROJ '/interim/BuoyancyFlux_c.mat']);
    flux = load(['../../Data/' PROJ '/processed/spurs1_flux_1hr.mat'],'Qr');
end
%ml_depth = ml_depth_04;

z_L = (dissipation.depth ./ L_pen')';
%z_L = (dissipation.depth ./ L)';

clrs = get( groot, 'defaultAxescolororder'); %COLORS
clrs(3,:) = clrs(4,:);

%TRY Interpolated Grid
t_mat = repmat(dissipation.time,[length(dissipation.depth),1]);
d_mat = repmat(dissipation.depth',[1,length(dissipation.time)]);

%interpolated grid
gi = ~isnan( log10(dissipation.epsilon));
[gridded_e] = griddata(t_mat(gi), d_mat(gi), log10(dissipation.epsilon(gi)), t_mat, d_mat,'natural');

contour_levels = [-12.5, [-8:0.5:-5.5]];
subplot(2,1,1)
%[h,c]=contourf(dissipation.time, dissipation.depth,log10(dissipation.epsilon));
[h,c]=contourf(dissipation.time, dissipation.depth,gridded_e,contour_levels);
set(c,'LineColor','none')
shading interp; axis ij;
%cmap = parula(5);
cmap = zebrajet(5);
colormap(cmap)
cb1=colorbar; 
datetick('x'); ylabel('Depth [m]');
ylabel(cb1,'log_{10}(\epsilon)')
hold on
%plot(dissipation.time, ml_depth_04,'-k','linewidth',2.5)
%plot(dissipation.time, ml_depth_04,'-w','linewidth',1)
plot(dissipation.time, movmean(ml_depth_04,24,'omitnan'),'-k','linewidth',2.5)
plot(dissipation.time, movmean(ml_depth_04,24,'omitnan'),'-w','linewidth',1)


%Plot zoom-box
if strcmp(PROJ,'SPURS-1')
    x1 = datenum('5/8/2013');
    x2 = datenum('5/24/2013');
    y1=0;
    y2 = 85;
elseif strcmp(PROJ,'SPURS-2')
    x1 = datenum('4/28/2017');
    x2 = datenum('5/10/2017');
    y1=0;
    y2 = 60;
end
plot([x1,x1],[y1,y2],'--r','linewidth',2)
plot([x1,x2],[y1,y1],'--r','linewidth',2)
plot([x2,x2],[y1,y2],'--r','linewidth',2)
plot([x1,x2],[y2,y2],'--r','linewidth',2)

%Plot Measurement locations:
plot(dissipation.time(1).*ones(length(dissipation.depth)),dissipation.depth,...
    'ko','linewidth',2,'markersize',8,'markerfacecolor',[1,1,1]*0.8)

title(PROJ)
caxis([-8,-5.5])

subplot(2,1,2)
contourf(t_mat, d_mat, gridded_e,contour_levels,'LineColor','none');
shading interp; axis ij;
cb1=colorbar; 
ylabel(cb1,'log_{10}(\epsilon)')
hold on
plot(dissipation.time, ml_depth_04,'-k','linewidth',2.5)
plot(dissipation.time, ml_depth_04,'-w','linewidth',1)
%Plot Measurement locations:
plot(datenum('5/8/2013').*ones(length(dissipation.depth)),dissipation.depth,...
    'ko','linewidth',2,'markersize',8,'markerfacecolor',[1,1,1]*0.8)

set(gca,'xlim',[x1,x2],'ylim',[y1,y2])
datetick('x','keeplimits'); ylabel('Depth [m]');
caxis([-8,-5.5])

%zoom-line annotations
annotation('line',[.59 .13],[.72 .45],'color','r','linewidth',2)
annotation('line',[.62 .8545],[.72 .45],'color','r','linewidth',2)

%subplot labels
annotation('Textbox',[0.132,0.9,0.015,0.02],'String','a',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
annotation('Textbox',[0.132,0.42,0.015,0.02],'String','b',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')

SPURS_figsave(gcf)

%% Wake Model (Appendix Figure)
% Simple monochromatic wake tests

clear
close all

a0 = 1;
g = 9.81;

f = 1/8; %freq. in hertz
omega = f*2*pi;
k = (omega).^2./g;

Ux = a0*omega*exp(-10*k)*0.075; %x times the wave orbitals
Uy = a0*omega*exp(-10*k)*0.075;

z = [-1:0.01:10]./k;
x = z;
y = x;
[XG,YG] = meshgrid(x,y);

z_aqd_initial = 10;

%initialize
pos_aqd = [0,0,z_aqd_initial]; %x, z
beam_vec =[0,1,0];

x_wake = 0;
y_wake = 0;
z_wake = z_aqd_initial+cos(-omega);

eta = [];
dt = omega^(-1) / 5; %ten dt per cycle

CamPos = [37.5887  -34.0101    2.3871];

figure(1),clf

beam_vec_past = pos_aqd;
for t = 0:dt:(4*(2*pi)/omega)
    sse = a0*cos(k.*XG - omega.*t);
    pos_aqd(3) = z_aqd_initial + cos(-omega.*t);
    az = a0*exp(-k.*z_wake);
    
    beam_vec = -[Ux+a0.*exp(-k.*pos_aqd(3)).*omega.*cos(-omega.*t),Uy,0];
    beam_vec = beam_vec./ sqrt(beam_vec(1)^2 + beam_vec(2)^2)*2 + pos_aqd;
    
    z_wake = cat(1, pos_aqd(3), z_wake + (az.*omega.*sin(k.*x_wake-omega.*t) )*dt);
    x_wake = cat(1, pos_aqd(1), x_wake + (Ux + az.*omega.*cos(k.*x_wake-omega.*t) )*dt);
    y_wake = cat(1, pos_aqd(2), y_wake + Uy*dt);
    
    inv_time = fliplr([-dt:dt:t]);
%     TEST 2-wave system
%     sse = a0*cos(k.*XG - omega.*t) + a2*cos(k.*XG - omega2.*t + phase);
%     pos_aqd(3) = z_aqd_initial + a0*cos(-omega.*t) + a2*cos(-omega2.*t + phase);
%     az = a0*exp(-k.*z_wake) + a2*exp(-k2.*z_wake);
%     
%     total_U = Ux+a0.*exp(-k.*pos_aqd(3)).*omega.*cos(-omega.*t) + a2.*exp(-k2.*pos_aqd(3)).*omega2.*cos(-omega2.*t + phase).*cos(theta2);
%     total_V = Uy + a2.*exp(-k2.*pos_aqd(3)).*omega2.*cos(-omega2.*t + phase).*sin(theta2);
%     
%     beam_vec = -[total_U,total_V,0];
%     beam_vec = beam_vec./ sqrt(beam_vec(1)^2 + beam_vec(2)^2)*2 + pos_aqd;
%     
%     az1 = a0*exp(-k.*z_wake);
%     az2 = a2*exp(-k2.*z_wake);
%     z_wake = cat(1, pos_aqd(3), z_wake + (az1.*omega.*sin(k.*x_wake-omega.*t) )*dt + az2.*(omega2.*sin(k2.*x_wake-omega2.*t+phase) )*dt);
%     x_wake = cat(1, pos_aqd(1), x_wake + (Ux + az1.*omega.*cos(k.*x_wake-omega.*t) )*dt + az2.*(omega2.*cos(k2.*x_wake-omega2.*t+phase) ).*cos(theta2)*dt);
%     y_wake = cat(1, pos_aqd(2), y_wake + Uy*dt + az2.*(omega2.*cos(k2.*y_wake-omega2.*t+phase) ).*sin(theta2)*dt);

    figure(1)
    subplot(2,1,1)
    surf(XG, YG, sse);
    shading flat
    set(gca,'zlim',[-3,3],'xlim',[-3,3],'ylim',[-3,3],'zdir','reverse',...
        'CameraPosition',[61.3119   -7.0080  -19.5487])
    
    subplot(2,1,2)
    plot3([0,0],[0,0],[a0*cos(-omega.*t),max(z)],'-k')
    hold on
    %plot3(x_wake, y_wake, z_wake,'-ok')
    scatter3(x_wake, y_wake, z_wake,10,inv_time,'filled')
    plot3(pos_aqd(1),pos_aqd(2),pos_aqd(3),'rx','linewidth',3)
    plot3([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)], [pos_aqd(3), beam_vec(3)],'r-','linewidth',3)
    hold off
    
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('Depth [m]')
    grid on
    set(gca,'zlim',[8,12],'xlim',[-3,3],'ylim',[-3,3],'zdir','reverse',...
        'CameraPosition',[61.3119   -7.0080  -19.5487])
    drawnow
    
    %
    if t > (3*(2*pi)/omega) && t < (3*(2*pi)/omega)+2*dt
        figure(15)
        subplot(2,4,1)
    
       plot3([0,0],[0,0],[a0*cos(-omega.*t),max(z)],'-k')
       hold on
       %plot3(x_wake, y_wake, z_wake,'-ok')
       scatter3(x_wake, y_wake, z_wake,10,inv_time,'filled')
       plot3(pos_aqd(1),pos_aqd(2),pos_aqd(3),'rx','linewidth',3)
       plot3([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)], [pos_aqd(3), beam_vec(3)],'r-','linewidth',3)
       hold off
       
       xlabel('x [m]')
       ylabel('y [m]')
       zlabel('Depth [m]')
       grid on
       set(gca,'zlim',[8,12],'xlim',[-3,3],'ylim',[-3,3],'zdir','reverse',...
           'CameraPosition',CamPos)
       drawnow
       
       subplot(2,4,5)
       hold on
       scatter(x_wake, y_wake,10,inv_time,'filled')
       plot(pos_aqd(1),pos_aqd(2),'rx','linewidth',3)
       plot([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)],'r-','linewidth',3)
       hold off
       xlabel('x [m]')
       ylabel('y [m]')
       grid on
       set(gca,'xlim',[-3,3],'ylim',[-3,3])
       
    elseif t > ((3+1/4)*(2*pi)/omega) && t < ((3+1/4)*(2*pi)/omega)+dt
        figure(15)
        subplot(2,4,2)
            
       plot3([0,0],[0,0],[a0*cos(-omega.*t),max(z)],'-k')
       hold on
       %plot3(x_wake, y_wake, z_wake,'-ok')
       scatter3(x_wake, y_wake, z_wake,10,inv_time,'filled')
       plot3(pos_aqd(1),pos_aqd(2),pos_aqd(3),'rx','linewidth',3)
       plot3([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)], [pos_aqd(3), beam_vec(3)],'r-','linewidth',3)
       hold off
       
       xlabel('x [m]')
       ylabel('y [m]')
       zlabel('Depth [m]')
       grid on
       set(gca,'zlim',[8,12],'xlim',[-3,3],'ylim',[-3,3],'zdir','reverse',...
           'CameraPosition',CamPos)
       drawnow
       
       subplot(2,4,6)
       hold on
       scatter(x_wake, y_wake,10,inv_time,'filled')
       plot(pos_aqd(1),pos_aqd(2),'rx','linewidth',3)
       plot([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)],'r-','linewidth',3)
       hold off
       xlabel('x [m]')
       ylabel('y [m]')
       grid on
       set(gca,'xlim',[-3,3],'ylim',[-3,3])
       
    elseif t > ((3+2/4)*(2*pi)/omega) && t < ((3+2/4)*(2*pi)/omega)+dt
        figure(15)
        subplot(2,4,3)
            
       plot3([0,0],[0,0],[a0*cos(-omega.*t),max(z)],'-k')
       hold on
       %plot3(x_wake, y_wake, z_wake,'-ok')
       scatter3(x_wake, y_wake, z_wake,10,inv_time,'filled')
       plot3(pos_aqd(1),pos_aqd(2),pos_aqd(3),'rx','linewidth',3)
       plot3([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)], [pos_aqd(3), beam_vec(3)],'r-','linewidth',3)
       hold off
       
       xlabel('x [m]')
       ylabel('y [m]')
       zlabel('Depth [m]')
       grid on
       set(gca,'zlim',[8,12],'xlim',[-3,3],'ylim',[-3,3],'zdir','reverse',...
           'CameraPosition',CamPos)
       drawnow
       
       subplot(2,4,7)
       hold on
       scatter(x_wake, y_wake,10,inv_time,'filled')
       plot(pos_aqd(1),pos_aqd(2),'rx','linewidth',3)
       plot([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)],'r-','linewidth',3)
       hold off
       xlabel('x [m]')
       ylabel('y [m]')
       grid on
       set(gca,'xlim',[-3,3],'ylim',[-3,3])
    elseif t > ((3+3/4)*(2*pi)/omega) && t < ((3+3/4)*(2*pi)/omega)+dt
        figure(15)
        subplot(2,4,4)
            
       plot3([0,0],[0,0],[a0*cos(-omega.*t),max(z)],'-k')
       hold on
       %plot3(x_wake, y_wake, z_wake,'-ok')
       scatter3(x_wake, y_wake, z_wake,10,inv_time,'filled')
       plot3(pos_aqd(1),pos_aqd(2),pos_aqd(3),'rx','linewidth',3)
       plot3([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)], [pos_aqd(3), beam_vec(3)],'r-','linewidth',3)
       hold off
       
       xlabel('x [m]')
       ylabel('y [m]')
       zlabel('Depth [m]')
       grid on
       set(gca,'zlim',[8,12],'xlim',[-3,3],'ylim',[-3,3],'zdir','reverse',...
           'CameraPosition',CamPos)
       drawnow
       
       subplot(2,4,8)
       hold on
       scatter(x_wake, y_wake,10,inv_time,'filled')
       plot(pos_aqd(1),pos_aqd(2),'rx','linewidth',3)
       plot([pos_aqd(1),beam_vec(1)], [pos_aqd(2), beam_vec(2)],'r-','linewidth',3)
       hold off
       xlabel('x [m]')
       ylabel('y [m]')
       grid on
       set(gca,'xlim',[-3,3],'ylim',[-3,3])
    end

    beam_vec_past = beam_vec;
end

% annotations---------------------------------------------

figure(15)
set(gcf,'Position',[34,85,1279,720])
x_left = 0.13;
y_top = 0.9;
y_bot = 0.42;

a = annotation('Textbox',[x_left, y_top,0.015,0.025],'String','a',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
b = annotation('Textbox',[0.34, y_top,0.015,0.025],'String','b',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
c = annotation('Textbox',[0.545, y_top,0.015,0.025],'String','c',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
d = annotation('Textbox',[0.75, y_top,0.015,0.025],'String','d',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')

e = annotation('Textbox',[x_left, y_bot,0.015,0.025],'String','e',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
f = annotation('Textbox',[0.34, y_bot,0.015,0.025],'String','f',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
g = annotation('Textbox',[0.545, y_bot,0.015,0.025],'String','g',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')
h = annotation('Textbox',[0.75, y_bot,0.015,0.025],'String','h',...
    'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',16,'BackgroundColor','w')

%% Call Shell Script to move new figures to Dropbox/Overleaf

!./mv_figs2.sh

%% FIGURE SAVE FUNCTION
function [flag] = SPURS_figsave(fighandle)
save_path = '../../Docs/Figures/PaperFigures2/';
figure_number_str =  num2str(get(fighandle,'Number'),'%02.f');

export_fig([save_path 'Figure' figure_number_str],'-pdf');

end