% Processing
%
% (1) (Check Raw)
%   if not ascii 2 mat
%
% (2) Check unwrap
%   if not call snaphu
%
% (2) Esimate Dissipations
%   (b)use wavenum Spectra ... need to simplify!
%
% (3) Save results to interim bursts
%
%
%

clear
tic

%FILE INFO
PROJ = 'SPURS-2';

if strcmp(PROJ,'SPURS-2')
    id_list = {'8089','6790','8116','5347'};%,'5143'}; %SPURS-2 (NOTE 5143 battery failed, no data return)
    %id_list = {'8089','6790','8116'};
    %id_list = {'6790'}; %SPURS-2
    %%NOW FIND FREQ WITH HEADER FILE
    %Freq_list = 2; 
    %Freq_list = [2,2,2,1,2]; %most are 2 MHz, except 5347
    
elseif strcmp(PROJ,'SPURS-1')
    id_list = {'5143','6774','8089','9134','6790','6524','8116'};
    %id_list = {'9134','6790','6524','8116'};  %didn't  finish first time...
    %id_list = {'6790'};
    %Freq_list = 2*ones(7,1); %NOW FIND WITH HEADER FILE
    %Freq_list(4) = 1;
    %name_list = {'c_000001','677404','808906','913410','679007','652404','811603'};
end

%PROCESSING OPTIONS
read_raw = 0; %1 to force re-read of raw files
force_unwrap = 0;%force unwrap
concat_bursts = 1; %after main processing loop, gather all dissipations into time series
plots = 0; %plot wavenumber spectra fits 0-off, 1-on

% LOOP THROUGH INSTRUMENT IDS
for HR_id = 1:length(id_list)
    
    Nortek_fn = id_list{HR_id};
    
    disp(['Instrument ID: ' Nortek_fn '...'])
    
    %GEN file pathing
    fp = ['../../Data/' PROJ '/interim/NortekFiles/' Nortek_fn '/'];
    raw_path = ['../../Data/' PROJ '/raw/NortekFiles/' Nortek_fn '/'];
    flag_path = ['../../Data/' PROJ '/interim/advection_flag.mat'];
    
    % NORTEK meta - check .hdr file to make sure these match w/ instrumnet
    % settings
    
    %THIS IS A PATCH HERE - I USED "instrument" while other code uses
    %"Head"
    try %SOME .hdr files are missing, fill in with default settings
        head_dir = dir([raw_path '*.hdr']);
        [Head] = readNortekHeader([raw_path head_dir(1).name]);
        
%         instrument.Vr = 0.06;
        instrument.nbins = Head.ncells;
        %Find xmit and rcv windows in units of length. Nortek defines in
        %vertical coords, so adjust to be along-beam as well with 25deg
        %angle
        if Head.freq == 2*10^6
            instrument.Lxmit = Head.Txmit*1500/2;
            instrument.Lrecv = Head.Trecv*1500/2;
            instrument.binSize = Head.Trecv*1500/2; %in meters
        elseif Head.freq == 1*10^6 %Clock cycles defined differently btwn 1 and 2 MHz instruments
            instrument.Lxmit = Head.Txmit*1500;
            instrument.Lrecv = Head.Trecv*1500;
            instrument.binSize = Head.Trecv*1500; %in meters
        end
        instrument.mAvg = Head.pings_per_sample; %number of avg'd pings
        instrument.range = [1:instrument.nbins]*instrument.binSize + Head.blank_dist;
        instrument.freq = Head.freq;
        instrument.sample_rate = Head.sample_rate;
    catch
        %instrument.Vr = 0.0632; %wrapping velocity noted in .hdr file (with added prec)
        instrument.range = [0.13:0.03:1.99]; %distance to bin  centers as  noted in .hdr file
        instrument.nbins = 63;
        instrument.freq = 2*10^6; %in hertz
        instrument.binSize = 0.027; %in meters
        instrument.Lxmit = 0.034;
        instrument.Lrecv = 0.027;
        instrument.mAvg = 10; %number of avg'd pings
        instrument.pulse_dist = 2.02;
        instrument.sample_rate = 8;
    end
    
    mean_c_cutoff = 60; %changed to 60 4/2/21 SFZ
    individual_c_cutoff = 40;
    qc_opts = 0; %0 - use rotations, 1 - don't use rotations
    [~,start_bin] = min( abs(instrument.range - 0.65) ); %automate at 30cm?
    %wavenumber index 1 is the zero frequency, index 2 seems to be biased
    %due to window effects?
    k_inds = 3:floor((length(instrument.range)-start_bin-1)/2 -4);
    taper = 1;
    noise_model = 'Fit'; %'Zedel1996'; %
    
    %Nortek_name = name_list{HR_id}; %Changed to scan for .v1 files
    %v1_files = dir([raw_path '*.v1']);
    v1_files = dir([raw_path '*.hdr']);
    for name_id = 1:length(v1_files)
        
        %fn = v1_files(name_id).name(1:end-3); %modified to -4 when using .hdr as reference file
        fn = v1_files(name_id).name(1:end-4);
        disp(['Instrument file: ' fn '...'])
        
        %% Check for burst files
        raw_files = dir([fp fn '_burst*.mat']);
        
        if isempty(raw_files)  || read_raw == 1
            disp('Converting Nortek Ascii files to .mat...')
            NortekHRAscii2burstMat(raw_path, fn, fp, instrument.nbins);
            disp('Finished converting raw files...')
            
            raw_files = dir([fp fn '_burst*.mat']);
        end
              
        %% Parfor loop
        parfor fi = 1:length(raw_files)
%        for fi  = 1:length(raw_files)
            %---------------- Load burst --------------------------------
            raw = load([fp raw_files(fi).name],'v1','a1','c1','sen','vUnwrap');
            C = struct(); %structure to hold processing results
            temp = struct();%structure to hold unsaved results
            
            %FIX MISTAKE I MADE SAVING EXTRA VARS, THIS TAKES UP TOO MUCH
            %SPACE... WILL BE SLOW
            %rmvar([fp raw_files(fi).name],'v1_demean','vQC','cQC')
            
            %display status
            if mod(fi,500)==0
                disp(['Burst ' num2str(fi) ' of ' num2str(length(raw_files))]);
            end
            
            %---------------- Unwrap --------------------------------
            %must have at least 10  pts to unwrap
            [len, nbins] = size(raw.v1);
            
            if (len >9 && ~isfield(raw,'vUnwrap')) || (force_unwrap == 1 && len >9)

                %March 8th, 2020 (SFZ) - replace SNAPHU w/ Una Millers histogram unwrap codes
                %[C.vUnwrap] = unwrap_UNA(raw.v1(:,3:end)');
                %C.vUnwrap = C.vUnwrap';
                
                %APRIL 20th -  SWITCH TO TOM's UNWRAP CODE
                C.Vr = max(max(abs(raw.v1(:,3:end))));
                temp.v1_demean  = raw.v1(:,3:end)-mean(raw.v1(:,3:end),2);
                C.vUnwrap = histogram_unwrap_function5(temp.v1_demean,C.Vr,2);
                
                %NOTE, raw files have 1st two columns as  burstnum, ensemblenum
                %hash = num2str(fi);
                %[C.vUnwrap] = snaphu2Dunwrap(raw.v1(:,3:end), raw.a1(:,3:end), raw.c1(:,3:end), instrument.Vr, pwd, hash);
                %C.vUnwrap = C.vUnwrap'; %Codes used to row/column swap...
                
                [len,nbins] = size(C.vUnwrap);
            elseif ~isfield(raw,'vUnwrap')
                C.vUnwrap = nan( size(raw.v1(:,3:end)) );
                C.Vr = nan;
                temp.v1_demean = C.vUnwrap;
            else
                C.vUnwrap = raw.vUnwrap;
            end
            
            %---------------- hash sensor file --------------------------------
            C.time = datenum(raw.sen(:,3), raw.sen(:,1), raw.sen(:,2), raw.sen(:,4), raw.sen(:,5), raw.sen(:,6));
            C.heading = raw.sen(:,13);
            %C.pitch =  raw.sen(:,14);
            %C.roll  = raw.sen(:,15);
            %C.pressure  = raw.sen(:,16);
            %C.temperature = raw.sen(:,17);
            
            %------------------------- Check to process dissipation --------------------------
            should_I_process = isfield(C,'vUnwrap') & len>9 & any(~isnan(C.vUnwrap));
            
            % call dissipation code
            if should_I_process
                
                %March 8th, 2020 (SFZ) - USE QC Function, rather than QC
                %inside wavenumber_dissipation
                
                %SOME .c1 files have leading burst num, some do not... use
                %adaptive indexing to solve
                %if both have leading column, ignore first 2 columns
                if size(raw.c1,2) == size(raw.v1,2)
                    C.col_inds = 3:nbins;
                elseif size(raw.c1,2) < size(raw.v1,2) %here corr does not have leading column
                    C.col_inds = 1:nbins;
                end
                
                %make sure same number of profiles
               if size(raw.c1,1) == size(raw.v1,1)
                    C.row_inds = 1:size(raw.c1,1);
                elseif size(raw.c1,1) < size(raw.v1,1)%here corr does not have leading column
                    C.row_inds = 1:nbins;
               elseif size(C.heading,1) < size(C.vUnwrap,1)
                   C.row_inds = 1:size(C.heading,1);
               end
                
                [temp.vQC, temp.cQC] = QC_nortek_velocities(C.vUnwrap, raw.c1(C.row_inds,C.col_inds),...
                    C.heading(C.row_inds), instrument.sample_rate, mean_c_cutoff, individual_c_cutoff, start_bin,qc_opts);
                
%                 if length(raw.c1(:,3:end)) == length(C.vUnwrap)
%                     [C.vQC, C.cQC] = QC_nortek_velocities(C.vUnwrap, raw.c1(:,end-nbins+1:end),[C.heading(1:len),C.pitch(1:len),C.roll(1:len)], mean_c_cutoff, individual_c_cutoff, start_bin,qc_opts);
%                 else
%                     [C.vQC, C.cQC] = QC_nortek_velocities(C.vUnwrap, raw.c1(1:length(C.vUnwrap),3:end),...
%                         [C.heading(1:length(C.vUnwrap)),C.pitch(1:length(C.vUnwrap)),C.roll(1:length(C.vUnwrap))], mean_c_cutoff, individual_c_cutoff, start_bin,qc_opts);
%                 end
%                 
                try
                [C.epsilon, C.epsilon_CI, C.Pxx, C.k, C.N, C.DOF, C.diss_flag, C.processing_opt_args] = wavenumber_dissipation(temp.vQC, temp.cQC, instrument, ...
                                                                            plots, taper, k_inds, noise_model);
                catch
                    C.epsilon = nan;
                    C.epsilon_CI = [nan,nan];
                    C.Pxx = nan;
                    C.k = nan;
                    C.N = [nan,nan];
                    C.DOF = nan;
                    C.diss_flag = nan;
                    C.processing_opt_args = nan;
                    %warning('wavenumber_dissipation failed...')
                end
            else
                C.epsilon = nan;
                C.epsilon_CI = [nan,nan];
                C.Pxx = nan;
                C.k = nan;
                C.N = [nan,nan];
                C.DOF = nan;
                C.diss_flag = nan;
                C.processing_opt_args = nan;
            end
            
            if C.diss_flag==2 %nan b/c low correlations
                C.flag = 2;
            else
                C.flag = 0; %should have exited normally and have a dissipation
            end
            %ADD META DATA LABEL .... maybe create new if appended??
            C.label = ['File creation time: ' datestr(now) ' Source Code: ' extractAfter(mfilename('fullpath'),'Projects') ];
            
            %---------------- append .mat files --------------------------------
            parsave([fp raw_files(fi).name],C);
        end
    end
    toc
    
    %% Concatenate bursts to make dissipation time series? Maybe do this in separate code??
disp('Concatenate bursts...')
tic
    if concat_bursts == 1
        raw_files = dir([fp '*_burst*.mat']);
        
        %find length of k
        ki=1;
        while ki ~= 0
            burst = load([fp raw_files(ki).name],'Pxx','k'); %need to allow variable length of Spectrum for 1/2MHz systems
            k_len = length(burst.k);
            if k_len==1
               ki=ki+1; 
            elseif k_len>1
                ki=0;
            end
            if ki>length(raw_files)
                disp('No Processed Spectra Found...')
                [fp raw_files(ki-1).name]
                k_len = 1;
                break
            end
            
        end
        diss.epsilon = nan(length(raw_files),1);
        diss.epsilon_CI = nan(length(raw_files),2);
        diss.time = nan(length(raw_files),1);
        diss.label = ['File creation time: ' datestr(now) ' Source Code: ' extractAfter(mfilename('fullpath'),'Projects') ];
        diss.k = nan(length(raw_files),k_len);
        diss.PSD = nan(length(raw_files),k_len);
        diss.N = nan(length(raw_files),2);
        diss.DOF = nan(length(raw_files),1);
        diss.instrument = instrument;
        
        %raw_files = dir([fp '*_burst*.mat']);
        for fi = 1:length(raw_files)
            %---------------- Load burst --------------------------------
            burst = load([fp raw_files(fi).name],'epsilon','epsilon_CI','time','Pxx','N','k','DOF');
            
            diss.epsilon(fi) = burst.epsilon;
            diss.epsilon_CI(fi,:) = burst.epsilon_CI;
            diss.time(fi) = mean(burst.time);
            if length(burst.Pxx)==1
                diss.PSD(fi,:) = nan(1,k_len);
                diss.k(fi,:) = nan(1,k_len);
            elseif length(burst.Pxx)>k_len
                %not sure why this would happen... seems like if a burst 
                %didn't get reprocessed.... maybe look into more later!!
                diss.PSD(fi,:) = nan(1,k_len);
                diss.k(fi,:) = nan(1,k_len);
                diss.epsilon(fi) = burst.epsilon;
            else
                diss.PSD(fi,:) = burst.Pxx;
                diss.k(fi,:) = burst.k;
            end
            diss.N(fi,:) = burst.N;
            diss.DOF(fi) = burst.DOF;
            
        end
        
        %sort in time
        [~, sort_ind] = sort(diss.time);
        diss.time = diss.time(sort_ind);
        diss.epsilon = diss.epsilon(sort_ind);
        diss.epsilon_CI = diss.epsilon_CI(sort_ind,:);
        diss.PSD = diss.PSD(sort_ind,:);
        diss.N = diss.N(sort_ind,:);
        diss.DOF = diss.DOF(sort_ind);
        diss.k = diss.k(sort_ind,:);
        
        save([ fp Nortek_fn '_dissipation_timeseries.mat'],'-struct','diss')
    end
    
toc
end

%GRID ALL - MAKE SURE  TO CHANGE  FILE  NAMES HERE!!!
grid_dissipation_data();