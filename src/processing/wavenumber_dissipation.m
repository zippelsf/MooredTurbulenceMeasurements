% CODE TO ESTIMATE THE TKE DISSIPATION RATE OF A NORTEK ADCP USING THE
% QUALITY CONTROLLED VELOCITY AND CORRELATIONS. THE CODE:
%   -CALCULATES THE POWER SPECTRUM,
%   -CORRECTS FOR THE NORTEK SPATIAL FILTER,
%   -OPTIONALLY ESTIMATES THE INSTRUMENT NOISE FROM CORRELATIONS
%   -DOES REGRESSIVE FIT TO A MODEL SPECTRUM
%   -RETURNS A TKE DISSIPATION RATE, NOISE ESTIMATE, AND THE POWER SPECTRUM
%
% INPUTS
%       vUnwrap         - [Nbins x nPts] double of velocity data
%       c1              - [nPts x Nbins+2] double (FIX THIS IN FUTURE, indexing is frustrating)
%       instrument      - Structure containing the fields listed below:
%           .pulse_dist - pulse distance
%           .freq       - frequency of instrument, used to estimate wrapping velocity
%           .mAvg       - Number of pings averaged, used in error estimate
%           .binSize    - in meters, used for wavenumber spacing
%           .Vr         - wrapping velocity in m/s
%       varargin        -
%           1  = ploting flag (1 = plot, other = noplot)
%           2  = Taper Options
%           3  = index of spectra to fit
%           4  = Noise estimation method (Fit, Zedel or Shcherbina)
%
% OUTPUTS
%       epsilon         - [1x1 double] TKE dissipation rate
%       epsilon_CI      - [1x2 double] TKE dissipation rate 95%CI estimates
%       P_avg           - [1xnfft double] Power Spectrum in (m^2 / k)
%       k               - [1xnfft double] Wavenumbers, k in [rad/m]
%       sigma2_k        - [1x1 double] White Noise estimate in  (m^2 / k)
%       DOF             - [1x1 double] Number of DOF take as 2*nPts
%       FLAG            - [1x1 double] Output flag (see below)
%       optargs         - TBD
%
% Code developed by S. Zippel. 2018/2019/2020. Notes:
%   -QC moved to separate function.
%   -Added options for regression (hardcoded)
%   -Added exptra outputs (spectrum, and noise estimates)
%
%FLAG reports if finished normally, or exited for QC reasons
%   1 - exit normally
%   2 - not enough good profiles (low correlations)
%  -9 - TOO MANY FUNCTION INPUTS (user error)                                                                         
                                                                            
function [epsilon, epsilon_CI, P_avg, k, sigma2_k, DOF, flag, optargs] = wavenumber_dissipation(vQC, cQC, instrument, varargin)

% set QC params
numvarargs = length(varargin);
if numvarargs > 4
    flag = -9;
    error('too many function inputs')
end

%1  = ploting flat (1 = plot, other = noplot)
%2  = taper flat (1 = hamming taper, 0 = no taper)
%3  = MIN_SNR, default = 4
%4  = Noise Model "Shcherbina2018" or "Zedel1996"

optargs = {1, 1, 4, 'Shcherbina2018'};
optargs(1:numvarargs) = varargin;

[plots, taper, k_inds, noise_model]  = optargs{:};
% plots = 1; %flag for plotting, 1=yes
% start_bin =  1;
% c_cutoff = 55;
% max_bad_pts_per_profile = 6; noise_model, 'Shcherbina2018')

epsilon = nan;
epsilon_CI = [NaN, NaN];
DOF = nan;
sigma2_k = nan;
k = nan;
P_avg = nan;
flag = 1; %FLAG is

[n_prfiles, nbins] = size(vQC);
DOF = 2*n_prfiles; %ROUGHLY

%% Calc Power spectrum for each good profile

nfft = nbins;
ks = 2*pi/instrument.binSize; %sampling wavenumber

if taper  == 1 %NOTE - MATLAB's periodogram seems to rescale for window gain already
    [pxx,k] = periodogram( detrend(vQC','linear',1), hamming(nfft),nfft, ks);
elseif taper == 0
    [pxx,k] = periodogram( detrend(vQC','linear',1), [],nfft, ks); %SWITCH TO USING LINEAR DETREND
end
P_avg  = mean(pxx,2);

%%  Estimate instrument noise from correlations

%NOISE PARAMS ESTIMATED FROM INSTRUMENT SETUP

if strcmp(noise_model, 'Shcherbina2018')
    M = instrument.mAvg; %number of avg'd samples
    sigma2_s  = ((cQC./100).^-2 - 1)./(2*M) .* ( instrument.Vr ./ pi).^2;%squared error earlier.. now same notation as in Shcherbina'18
    sigma2_s(isinf(sigma2_s))  =  nan; %in case  of R2=0
    sigma2_s_avg  =  nanmean(nanmean(sigma2_s));
    sigma2_k = sigma2_s_avg./k(end);

elseif strcmp(noise_model, 'Zedel1996')
    M = instrument.mAvg; %number of avg'd samples
    sigma2_s  = log((cQC./100)) .* (-2) .* ( instrument.Vr ./ pi).^2 ./M;
    sigma2_s(isinf(sigma2_s))  =  nan; %
    sigma2_s_avg  =  nanmean(nanmean(sigma2_s));
    sigma2_k = sigma2_s_avg./k(end);
end

%% Estimate dissipation

good_pts = zeros(size(k));
good_pts(k_inds) = 1;
good_pts = boolean(good_pts);

[epsilon, sigma2_k] = fit_ISR(P_avg(good_pts), k(good_pts),'Linear','LS', -99, instrument.Lxmit, instrument.Lrecv);
bad_pts = [];

%confidence estimates based loosely on chi2 (not tested)
sum_DOF = DOF*sum(good_pts);
epsilon_CI = epsilon .* ([chi2inv(1-0.95, sum_DOF)/sum_DOF, chi2inv(0.95, sum_DOF)/sum_DOF]).^(3/2);

% REMOVED - was erroniously left in after testing, 5/26/21
% %% QC Check (if fit SNR is larger than 95%CI)
% %coef1  = 18/55 * (8/9/0.4)^(2/3); %Veron&Melville 1999 JTECH coeff
% coef1 = 0.53; %changed to Sreenivasan 1995
% fit_snr = coef1*k(3).^(-5/3).*epsilon.^(2/3) ./ sigma2_k(1);
% 
% alpha = 0.01;
% err_high = DOF*sum(good_pts)./chi2inv(1-alpha/2, DOF*sum(good_pts));
% if (fit_snr < (1-err_high) )
%     epsilon = NaN;
% end

%% PLOTS -
if plots == 1
    % RESIDUALS/QC/ERROR ANALYSIS
    %coef1  = 18/55 * (8/9/0.4)^(2/3); %Veron&Melville1999
    chi_plus = chi2inv(0.95, DOF)/DOF;
    chi_minus = chi2inv(1-0.95, DOF)/DOF;

    figure,clf
    hold on
    plot(k(2:end), P_avg(2:end),'linewidth',2)
    plot(k(good_pts(bad_pts)), P_avg(good_pts(bad_pts)),'ro','linewidth',2,'HandleVisibility','off')
    plot(k, k.^(-5/3) * coef1 .* epsilon.^(2/3)  + sigma2_k,  '-k','linewidth',2)
    plot(k, (k.^(-5/3) * coef1 .* epsilon.^(2/3)  + sigma2_k)*chi_minus,  '-k','linewidth',1,'HandleVisibility','off')
    plot(k, (k.^(-5/3) * coef1 .* epsilon.^(2/3)  + sigma2_k)*chi_plus,  '-k','linewidth',1,'HandleVisibility','off')
    plot(k, ones(length(k),1)*sigma2_k,  '--r','linewidth',2)
    
    set(gca,'yscale','log','xscale','log')
    legend({'Spectral Est.',...
        ['\epsilon = ' num2str(epsilon,2)],'Noise Estimate'},...
        'fontsize',14,'location','best')
    grid  on
    xlabel('Wavenumber [rad/m]')
    ylabel('Power Spectral Density [m^2s^{-2} / k]')
end
end

%% Fitting fuction Added 2/25/2020 - try robustfit, have options under single
%fitting function
function [epsilon, Noise_est] = fit_ISR(Pxx, k, varargin)

%Options (linear or log sum of squares, robustfit)
%Specify Noise, or fit noise
% set QC params
numvarargs = length(varargin);
if numvarargs >5
    flag = -9;
    error('too many function inputs')
end

%1  = Fit Type (linear, log)
%2  = Fit Method (Least squares, robustfit)
%3  = Fit Noise or use est'd Noise
optargs = {'Linear', 'LS', 0, 0.036, 0.034};
optargs(1:numvarargs) = varargin;

[FitType, FitMethod, Noise_est, L1, L2]  = optargs{:};

if size(Pxx) ~= size(k)
    k = k';
end

%coef1  = 18/55 * (8/9/0.4)^(2/3); %Veron&Melville 1999 JTECH coeff
coef1  = 0.53; %changed to Sreenivasan et al 1995
     
if strcmp(FitType,'Linear') && strcmp(FitMethod,'LS') && (Noise_est == -99) %FIT PINK NOISE (??)
    %Fit with double noise (correlation and electronic)
    %Spectral model is P_meas = G*(Coef*e^(2/3)k^(-5/3) + N1)
    % G is gain function, N is noise due to correlations

    %Moved back to ML divide after Reviewer1 comments ... (SFZ 3/29/21)
    G = sinc(L1/2/pi.*k).^2 .* sinc(L2/2/pi.*k).^2;
    Y = Pxx;
    X = cat(2, G.*ones(size(Pxx)), G.*coef1.*k.^(-5/3));
    b = X\Y;
    model_fits(1) = b(2).^(3/2);
    model_fits(2) = b(1);
    
elseif strcmp(FitType,'Linear') && strcmp(FitMethod,'LS') && (Noise_est ~=0)
    %minimize Linear sum of squares w/ est'd noise
    min_func1 =@(model) sum( ((10^7)*Pxx - (10^7)*(Noise_est + coef1*abs(model(1)).^(2/3).*k.^(-5/3)) ).^2 );
    options = optimset('TolX',10^-7,'TolFun',10^-7,'FunValCheck','on','MaxIter',200);
    [model_fits, ~, exitflag] = fminsearch(min_func1, 10^-3,options);
    
elseif strcmp(FitType,'Linear') && strcmp(FitMethod,'LS') && (Noise_est ==0)
    min_func1 =@(model) sum( (Pxx - (abs(model(2)) + coef1*abs(model(1)).^(2/3).*k.^(-5/3)) ).^2 );
    options = optimset('TolX',10^-9,'TolFun',10^-9,'FunValCheck','on','MaxIter',250);
    [model_fits, ~, exitflag] = fminsearch(min_func1, [10^-3,10^-7],options);
    model_fits(2) = abs(model_fits(2));

elseif strcmp(FitType,'Log') && strcmp(FitMethod,'LS') && (Noise_est ~=0)
    %minimize Log sum of squares w/ est'd noise
    min_func1 =@(model) sum( (log(Pxx) - log(Noise_est + coef1*abs(model(1)).^(2/3).*k.^(-5/3)) ).^2 );
    options = optimset('TolX',10^-7,'TolFun',10^-7,'FunValCheck','on','MaxIter',200);
    [model_fits, ~, exitflag] = fminsearch(min_func1, 10^-3,options);

elseif strcmp(FitType,'Log') && strcmp(FitMethod,'LS') && (Noise_est ==0)
    %minimize Log sum of squares w/ fit noise
    min_func1 =@(model) sum( (log(Pxx) - log(model(2) + coef1*abs(model(1)).^(2/3).*k.^(-5/3)) ).^2 );
    options = optimset('TolX',10^-7,'TolFun',10^-7,'FunValCheck','on','MaxIter',200);
    [model_fits, ~, exitflag] = fminsearch(min_func1, [10^-3,10^-8],options);
    
elseif strcmp(FitType,'Linear') && strcmp(FitMethod,'Robust') && (Noise_est ==0)
    [b,stats] = robustfit( k.^(-5/3), Pxx);
    model_fits(2) = b(1);
    model_fits(1) = (b(2)./coef1).^(3/2);
    
elseif strcmp(FitType,'Linear') && strcmp(FitMethod,'Robust') && (Noise_est ~=0)
    [b,stats] = robustfit( k.^(-5/3), Pxx-Noise_est,'bisquare',4.685,'off');
    model_fits(1) = (b(1)./coef1).^(3/2);
end

epsilon = abs(model_fits(1));
if Noise_est ==0 || Noise_est == -99
    Noise_est = abs(model_fits(2:end));
end

end