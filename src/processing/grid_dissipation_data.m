% Script to make dissipation grid on the same time series as the flux
% hourly grid... use closest points
%
%
% REMEMBER TO UPDATE SAVE NAME

clear

PROJ = 'SPURS-2';
fp = ['../../Data/' PROJ '/interim/NortekFiles/'];

f_list = dir([fp '**/*timeseries.mat']);
if strcmp(PROJ,'SPURS-2')
    flux = load(['../../Data/' PROJ '/processed/spurs2_flux_1hr.mat'],'QN','mday','ustar');
    met  = load(['../../Data/' PROJ '/processed/spurs2_met_1hr.mat'],'prate','mday','sst');
    
    id_depth_order = {'8089','6790','8116','5347','5143',};
    aqd_depth = [7, 21.5, 41.5, 61.7, 101]; %in meters
    save_name = [fp 'SPURS2_dissipation_grid_v1c.mat'];
elseif strcmp(PROJ,'SPURS-1')
    flux = load(['../../Data/' PROJ '/processed/spurs1_flux_1hr.mat'],'QN','mday','ustar');
    
    id_depth_order = {'5143','6774','8089','9134','6790','6524','8116'};
    aqd_depth = [12.5, 21.5, 41.5, 61.7, 82, 101.6, 121.6]; %in meters
    %id_depth_order = {'5143','6774','8089'};
    %aqd_depth = [12.5, 21.5, 41.5,]; %in meters
    save_name = [fp 'SPURS1_dissipation_grid_v1c.mat'];
end

%PREALLOCATE
dissipation.time = flux.mday;
dissipation.depth = aqd_depth;
dissipation.epsilon = nan(length(dissipation.depth), length(dissipation.time));

%wrangle time series
for fi = 1:length(f_list)
    
    f_path  = [f_list(fi).folder '/' f_list(fi).name];
    
    diss_raw(fi) = load(f_path,'-mat');
    diss_raw(fi).epsilon(diss_raw(fi).epsilon < 10^-10) = NaN;
    diss_raw(fi).epsilon(diss_raw(fi).epsilon > 10^-4) = NaN;
    
    %Get water depth
    Nortek_ID = f_list(fi).folder(end-3:end);
    depth_index = find(strcmp(Nortek_ID, id_depth_order));
    
    %diss(fi).depth = aqd_depth(depth_index);
    
    %Get local flux data, calc z/L
    inds = dsearchn(diss_raw(fi).time, flux.mday);
    
    dissipation.epsilon(depth_index,:) = diss_raw(fi).epsilon(inds);
    dissipation.epsilon_95ci_lo(depth_index,:) = diss_raw(fi).epsilon_CI(inds,1);
    dissipation.epsilon_95ci_hi(depth_index,:) = diss_raw(fi).epsilon_CI(inds,2);
    dissipation.Spectra(depth_index).PSD = diss_raw(fi).PSD(inds,:);
    dissipation.Spectra(depth_index).k = diss_raw(fi).k(inds,:);
    dissipation.Spectra(depth_index).DOF = diss_raw(fi).DOF(inds);
    dissipation.Spectra(depth_index).N = diss_raw(fi).N(inds,:);
    
    %CHECK THAT TIMES ARE WITHIN AN HOUR OF EACH OTHER  
     t_diff = abs(flux.mday - diss_raw(fi).time(inds));
     bad_pts = t_diff > (0.5/24); %if time mismatch is greater than an hour, don't use...
     dissipation.epsilon(depth_index,bad_pts) = NaN;
     dissipation.epsilon_95ci_lo(depth_index,bad_pts) = NaN;
     dissipation.epsilon_95ci_hi(depth_index,bad_pts) = NaN;
     
     dissipation.Spectra(depth_index).PSD(bad_pts,:) = NaN;
     dissipation.Spectra(depth_index).k(bad_pts,:) = NaN;
     dissipation.Spectra(depth_index).DOF(bad_pts) = NaN;
     dissipation.Spectra(depth_index).N(bad_pts,:) =NaN;
     
     dissipation.meta{fi} = diss_raw(fi).label;
end


%% Save to /interim folder
[rows,cols] = size(dissipation.epsilon);
[t_rows,t_cols] = size(dissipation.time);

if t_cols ~= cols %match column dimension with time vector
    dissipation.time = dissipation.time';
end

save(save_name,'dissipation')
