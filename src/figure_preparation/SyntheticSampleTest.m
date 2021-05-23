% Figure 2 - synthetic sampling data
%
%
%

clear

N = 64;
oversample= 10;
n_profiles = 1000;
beta = -5/3;
clrs = colororder;

%design sample fitler
L1 = 12;%0.027/0.034*oversample; %Normalized window transmit and receive lengths
L2 = 0.034/0.034*oversample;

u = ones(floor(L1),1);
v = ones(floor(L2),1);
w = conv(u,v); %make convolution window
w = w./sum(w); %normalize window
filter1 = zeros(1,N*oversample);
filter1(1:length(w)/2) = w(length(w)/2+1:end);
filter1(end-length(w)/2+1:end) = w(1:length(w)/2);
filter1=filter1./sum(filter1);
len1 = length(u);

% % %filter of second length
% % L1 = 0.025/0.034*oversample;
% % u = ones(floor(L1),1);
% % w2 = conv(u,v); %make convolution window
% % w2 = w2./sum(w2); %normalize window
% % filter2 = zeros(1,N*oversample);
% % filter2(1:length(w2)/2) = w(length(w2)/2+1:end);
% % filter2(end-length(w2)/2+1:end) = w2(1:length(w2)/2);
% % len2 = length(u);
% % 
% % %Make time domain filter, Nx(N*oversample);
% % G = zeros(N,(N*oversample));
% % G2 = G; %G2 is the jittery sample filter
% % jitter_repeat = 10;
% % n1 = -1;
% % n2 = -1;
% % for ni = 1:N
% %    G(ni,:) = circshift(filter1,(ni-1)*len1);
% %    if mod(ni,jitter_repeat)==0
% %        n2 = n2+1;
% %        G2(ni,:) = circshift(filter2,n1*len1+n2*len2);
% %    else
% %        n1=n1+1;
% %        G2(ni,:) = circshift(filter1,n1*len1+n2*len2);
% %    end
% % end

% Create time series
yt = zeros(N*oversample,n_profiles);
yt_subfilt = zeros(N,n_profiles);
yt_jitter = yt_subfilt;
yt_jit_noise = yt_jitter;

for ii = 1:n_profiles
    [yt(:,ii)]= (synthetic_timeseries(N*oversample,beta) );
    yt2sample(:,ii) = yt(:,ii) + randn(size(yt(:,ii)))/2 ;
    %yt_subfilt(:,ii) = G*(yt2sample(:,ii));
    yt_conv(:,ii) = conv(yt2sample(:,ii),w,'same');
    
    %yt_subfilt(:,ii) = G*yt(:,ii);
%     yt_jitter(:,ii) = G2*yt(:,ii);  
end

yt_subsample = yt2sample(1:oversample:(N*oversample),:);
yt_subfilt = yt_conv(1:oversample:(N*oversample),:);
%yt_subsample = yt(1:oversample:(N*oversample),:);


%Periodograms
[Pyy,k1] = periodogram(yt, [],N*oversample,oversample);
[Pyy_sub,k2] = periodogram(yt_subsample, [],N,1);
[Pyy_subfilt,k2] = periodogram(yt_subfilt, [],N,1);
[Pyy_jitter,k2] = periodogram(yt_jitter, [],N,1);
% [Pyy_filter_func,k2] = periodogram(G, [],N,1);
[P_test] = periodogram(yt_conv, [], N*oversample, oversample);

Pyy = mean(Pyy,2);
Pyy_sub = mean(Pyy_sub,2);
Pyy_subfilt = mean(Pyy_subfilt,2);
% Pyy_jitter = mean(Pyy_jitter,2);
% Pyy_filter_func = mean(Pyy_filter_func,2);
P_test = mean(P_test,2);

%Plots
figure(2),clf
hold on
grid on
plot(k1(1:end-1)*N/2,Pyy(1:end-1),'linewidth',2)
plot([k2(2),k2(end)]*N/2, 0.25./(k2(end))*[1,1],'--k','linewidth',2)
plot(k2(1:end-1)*N/2,Pyy_sub(1:end-1),'linewidth',2,'color',clrs(2,:))
plot(k2(1:end-1)*N/2,Pyy_subfilt(1:end-1),'linewidth',2,'color',clrs(3,:))

%corrected plot
% L1 = 0.027/0.027*oversample; %Normalized window transmit and receive lengths
% L2 = 0.034/0.027*oversample;
response_func = sinc( L1/oversample*k2(1:end-1) ).^2 .* sinc( L2/oversample*k2(1:end-1) ).^2;
plot(k2(1:end-1)*N/2, Pyy_subfilt(1:end-1)./response_func,'linewidth',2,'color',clrs(4,:))
plot(k2(1:end-1)*N/2, response_func/1e2,'--','linewidth',2,'color',clrs(5,:))

%The below code will plot the full (not downsampled) effects of convolution
% response_func = sinc( L1/oversample*k1(1:end-1) ).^2 .* sinc( L2/oversample*k1(1:end-1) ).^2;
% plot(k1(1:end-1)*N/2, P_test(1:end-1)./response_func,'linewidth',2)
% plot(k1(1:end-1)*N/2,P_test(1:end-1),'linewidth',2)
% plot(k1(1:end-1)*N/2, P_test(1:end-1)./Pyy(1:end-1)*1e-2 )

%filter_power_spec = Pyy_filter_func ./ Pyy_filter_func(2); %normalize
%plot(k2(1:end-1)*N/2, Pyy_subfilt(1:end-1)./filter_power_spec(1:end-1),'linewidth',2)

set(gca,'yscale','log','xscale','log','ylim',[1e-3,1e2],'xlim',[0.3,3*1e2])
legend({'Original Signal','Added Noise','Uniform Subsample of Noisy Signal',...
    'Synthetic ADCP Sampling of Noisy Signal','Synthetic Sampling, Corrected','Response Function (offset by 10^{-2})'},'location','best')
xlabel('Wavenumber [cpm]')
ylabel('Velocity PSD [m^2 s^{-2} cpm^{-1}]')
