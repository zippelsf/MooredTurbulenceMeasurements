% Nortek line-by-line Read-in
%
%
%


function [] = NortekHRAscii2burstMat(fp, fn, savePath, nbins)

if nargin == 0 
    fp = '../../Data/SPURS-1/raw/NortekFiles/5143/';
    fn = 'c_000001';
    
    %nbins HR 2MHz AQDs  = 63, HR 1MHz AQDs = 70
    nbins = 63;
    
    %save info
    savePath = '../../Data/SPURS-1/interim/NortekFiles/5143/';
end

%Number of burst digits, here assume order 1000 (less than 9999 total, or 416 days)
burst_digits = 4;
name_zeros = ['%0' num2str(burst_digits) 'd'];

%Number of bins, here set to 63
beam_formatspec = repmat(['%f  '],[1, nbins+2]);
sen_formatspec = repmat(['%f  '],[1, 19]);

%Open files
fid_v1 = fopen([fp fn '.v1']);
fid_a1 = fopen([fp fn '.a1']);
fid_c1 = fopen([fp fn '.c1']);
fid_sen = fopen([fp fn '.sen']);

ii = 0; %Burst Counter
ei = 1;

%initialize vars
temp.v1 = [];
temp.a1 = [];
temp.c1  = [];
temp.sen = [];

while ~feof(fid_v1)
    
    ii  = ii+1; %while  loop counter
    while ei == ii %stop looping when data burst num changes
        
        V = textscan(fid_v1, beam_formatspec, 1);
        A = textscan(fid_a1, beam_formatspec, 1);
        C = textscan(fid_c1, beam_formatspec, 1);
        SEN = textscan(fid_sen, sen_formatspec, 1);
        
        ei = SEN{7}; %Grab Burst index from data ... MAKE SURE THE 7TH COLUMN IS BURSTNUM IN .HDR FILE
        
        %append burst data
        temp.v1 = cat(1,temp.v1,cell2mat(V));
        temp.a1 = cat(1,temp.a1,cell2mat(A));
        temp.c1 = cat(1,temp.c1,cell2mat(C));
        temp.sen = cat(1,temp.sen,cell2mat(SEN));
    end
    
    %prep vars to save as burst
    raw.v1 = temp.v1(1:end-1,:);
    raw.a1 = temp.a1(1:end-1,:);
    raw.c1 = temp.c1(1:end-1,:);
    raw.sen = temp.sen(1:end-1,:);
    
    %reset temp struct for new burst
    temp.v1(1:end-1,:) = [];
    temp.a1(1:end-1,:) = [];
    temp.c1(1:end-1,:) = [];
    temp.sen(1:end-1,:) = [];
    
    saveName = sprintf([savePath fn '_burst' name_zeros '.mat'],ii);
    save(saveName,'-struct','raw')
    
end
%     %TEST READIN! BURSTNUM (burstcounter) AND (ensembleCounter) INDEX-IN-BURST
%     check_burstnum_v1 = all(raw.v1(:,1) == ii);
%     check_burstnum_a1 = all(raw.a1(:,1) == ii);
%     check_burstnum_c1 = all(raw.c1(:,1) == ii);
%     check_burstnum_sen  = all(raw.sen(:,7)  == ii);
%     
%     check_ensemble_v1 = all(raw.v1(:,2) == [1:burst_len]');
%     check_ensemble_a1 = all(raw.a1(:,2) == [1:burst_len]');
%     check_ensemble_c1 = all(raw.c1(:,2) == [1:burst_len]');
%     check_ensemble_sen = all(raw.sen(:,8) == [1:burst_len]');
%     error_flag = any(~check_burstnum_v1 || ~check_burstnum_a1 || ~check_burstnum_c1 ||...
%         ~check_burstnum_sen || ~check_ensemble_v1 || ~check_ensemble_a1 ||  ~check_ensemble_c1 ||...
%         ~check_ensemble_sen);
%     
%     if error_flag
%         disp('Unexpected burstCounter or ensembleCounter...')
%         return
%     end

fclose(fid_v1);
fclose(fid_a1);
fclose(fid_c1);
fclose(fid_sen);

% 
%     if feof(fid_v1)
%         fclose(fid_v1);
%         fclose(fid_a1);
%         fclose(fid_c1);
%         fclose(fid_sen);
%         %disp('End of files...')
%         return
%     end
