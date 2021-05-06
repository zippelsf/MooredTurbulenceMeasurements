% Data QC for wavenumber processing
%
% (A) Remove first few bins for wake issues
% (B) Remove profiles with low mean correlations (mean_c_cutoff)
%       AND variance greater than 0.015
%       AND heading rate of change > 0.1821
% (C) Remove individual meas. with low correlations (individual_c_cutoff)
%       (AND) individual meas. larger than 3.75*stdev
% (D) Remove profiles with 3 or more consecutive NaNs
% (E) Interpolate over missing data linearly

% Inputs
%   v1                      [n_profiles x nbins] double of velocities
%   c1                      [n_profiles x nbins] double of correlations
%   heading                 [n_profiles x 3] double of instrument heading
%                           pitch and roll (in degs)
%   mean_c_cutoff           [1 x 1] double of mean profile corr cutoff
%   individual_c_cutoff     [1 x 1] double of individual meas corr cutoff
%   Opt_args                TBD
%
% Outputs
%   QC velocity matrix
%   correlation matrix

function [vQC, cQC, interp_percentage] = QC_nortek_velocities(v1, c1, heading, sample_rate, mean_c_cutoff, individual_c_cutoff, start_bin, Opt_args)

%Set Default Args
if nargin ==4  
    sample_rate = 8;
    mean_c_cutoff = 40;
    individual_c_cutoff = 40;
    start_bin = 6;
    Opt_args = 0;
elseif nargin == 5
    individual_c_cutoff = 45;
    start_bin = 6;
    Opt_args = 1;
elseif nargin == 6
    start_bin = 6;
    Opt_args = 1;
elseif nargin == 7
    Opt_args = 1;
end

%preallocate QC vars
vQC = NaN;
cQC = NaN;
interp_percentage = NaN;

%sample_rate = 8; %in Hertz SFZ-4/21/2020 NOW AN INPUT
var_cutoff = 0.015; %in (m/s)^2
%rotational_max = 4;%0.1821; %in rad/s should be roughly (sample_rate/2*binSZ/Rmax)
rotational_max = 1; %modified 3/29/21 SFZ

% if nargin ==3 
%     rotational_max = 0.4;
% end
%% STEP A - remove first n number of bins

v1 = v1(:,start_bin:end);
c1 = c1(:,start_bin:end);

%% STEP B - remove profiles with low mean c, high var, and low rotation

heading_rate_of_change = [(heading(2,1)-heading(1,1)); (heading(3:end,1) - heading(1:end-2,1))/2; (heading(end,1)-heading(end-1,1))].*sample_rate;
heading_rate_of_change = deg2rad(heading_rate_of_change);

mean_alongbeam_c = nanmean(c1,2);
mean_profile_var = var(v1,[],2);

if Opt_args == 0 %USE ALL SUGGESTED CUTOFFS
    good_profiles = find( mean_alongbeam_c > mean_c_cutoff &...
                        mean_profile_var < var_cutoff & ...
                        abs(heading_rate_of_change) < rotational_max);
elseif Opt_args == 1 %USE ONLY mean-cor and mean-var CUTOFFS
    good_profiles = find( mean_alongbeam_c > mean_c_cutoff &...
                        mean_profile_var < var_cutoff );
elseif Opt_args == 2 %USE ONLY mean-cor CUTOFF
    good_profiles = find( mean_alongbeam_c > mean_c_cutoff); 
end

if ~isempty(good_profiles)
    vQC = v1(good_profiles,:);
    cQC = c1(good_profiles,:);
else %NO GOOD PROFILES FOUND, END IMMEDIATELY
    disp('No Good Profiles found - mean profile correlations too low')
    return
end


%% Remove profiles with VARIANCE SPIKES

% mean_var = mean(var(vQC,[],2));
% var_inliers = find(var(vQC,[],2) < 3.5*mean_var);
% vQC = vQC(var_inliers,:);
% cQC = cQC(var_inliers,:);

%increased here, seems to reduce too much?
mean_var = mean(var(vQC,[],2));
var_inliers = find(var(vQC,[],2) < 3.5*mean_var);
vQC = vQC(var_inliers,:);
cQC = cQC(var_inliers,:);

% %repeat twice...
% mean_var = mean(var(vQC,[],2));
% var_inliers = find(var(vQC,[],2) < 3.5*mean_var);
% vQC = vQC(var_inliers,:);
% cQC = cQC(var_inliers,:);

%% STEP C - remove individual meas. w/ low c or spikes

vQC( cQC < individual_c_cutoff) = NaN;
cQC( cQC < individual_c_cutoff) = NaN;

% %REMOVE "SPIKES"
along_beam_stdev = nanstd(vQC,[],1);
%spikes = abs(vQC-nanmean(vQC,1) - 3.75*along_beam_stdev) > 1;
spikes = abs(vQC-nanmean(vQC,1)) > 3.75*along_beam_stdev;

vQC(spikes) = NaN;
cQC(spikes) = NaN;

%% STEP D - remove profiles with 3 or more consecutive NaNs, or less than 50% data return

[nProfiles, nBins] = size(vQC);
no_consec_nan_profiles = ones(1,nProfiles);
for ii = 1:nProfiles
    %search for sets of 3 consecutive NaNs. Returns empty if non found
    consec_nan = strfind( isnan(vQC(ii,:)), true(1,3) );
    if ~isempty(consec_nan) %|| sum(isnan(vQC(ii,:))) > (length(vQC(ii,:))/5) %found string of three NaNs, flag profile as bad
        no_consec_nan_profiles(ii) = 0;
    end
end

vQC = vQC(boolean(no_consec_nan_profiles),:);
cQC = cQC(boolean(no_consec_nan_profiles),:);

%% STEP E - Linear interp over NaN values

[nProfiles, nBins] = size(vQC);
x = 1:nBins;

interp_percentage = sum(sum(isnan(vQC)))./numel(vQC);

for ii=1:nProfiles
    temp_vel = vQC(ii,:);
    gi = ~isnan(temp_vel);
    
    %try replacing with mean value?
    %vQC(ii,~gi) = nanmean(vQC(ii,:));
    
    %TRY FILLGAPS - AR MODEL
    %vQC(ii,:) = fillgaps(vQC(ii,:),[],'aic');
    
    %vQC(ii,:) = interp1( x(gi), temp_vel(gi), x ,'nearest', 0);
    vQC(ii,:) = interp1( x(gi), temp_vel(gi), x ,'linear', 0);
    %vQC(ii,:) = interp1( x(gi), temp_vel(gi), x ,'pchip', 0);
end


end