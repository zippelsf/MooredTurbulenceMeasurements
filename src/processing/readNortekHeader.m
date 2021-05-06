function [Head] = readNortekHeader(infile)
%function [AquaDopp] = readNortekHeader(infile)


%  Reads in header information for Nortek AquaDopp HRP
%  
%  Input:  string containing file name (location + name without file extension)
%
%  Output:  Structure Head with pertinent header information 
%
%  Deborah Le Bel
%  02.05.2009
%  _____________________________________________________________________________



if nargin < 1
   help readNortekHeader;
   return;
end


initial_col = 39;


file = [infile];
fid = fopen(file);
for i=1:3
   dum = fgetl(fid);
end
nmeas = str2num(dum(initial_col:end));



for i=1:7
   dum = fgetl(fid);
end
for i=30:length(dum)-2
   if strcmp(dum(i:i+2),'sec')
      pos = i;
   end
end
burst_interval = str2num(dum(initial_col:pos-1));

dum = fgetl(fid)
for i=30:length(dum)-1
   if strcmp(dum(i:i+1),'m')	% changed to mm 02.04.2014
%   if strcmp(dum(i:i+1),'mm')
      pos = i;
   end
end
bin_size = str2num(dum(initial_col:pos-3));

if strcmp(dum(i:i+1),'mm')
   bin_size = bin_size/1000;
end



for i=1:9
   dum = fgetl(fid);
end
ncells = str2num(dum(initial_col:end)); 



for i=1:2
   dum = fgetl(fid);
end

for i=30:length(dum)
   if strcmp(dum(i),'m')
      pos = i;
   end
end
blank_dist =  str2num(dum(initial_col:pos-1));



for i=1:3
   dum = fgetl(fid);
end
samples_per_burst = str2num(dum(initial_col:end));



dum = fgetl(fid);
for i=1:length(dum)-1
   if strcmp(dum(i:i+1),'Hz')
      pos = i;
   end
end
sample_rate = str2num(dum(initial_col:pos-1));


for i=1:10
   dum = fgetl(fid);
end
nbeams = str2num(dum(initial_col:end));

dum = fgetl(fid);
pings_per_burst = str2num(dum(initial_col:end));

dum = fgetl(fid);
sver = str2num(dum(initial_col:42));

%Looks like some hardcoded correction for vertical vs. along-beam spacing?
%Not clear to me. I'll set the transmit and receive windows calculated as
%shown by Sven:
for i=1:5
   dum = fgetl(fid);
end
system1 = str2num(dum(initial_col:end));
Txmit = system1/125e3; %in us
%FROM SVEN
% The durations here might be somewhat higher than you would expect. 
%This is because the AquaProHR software always calculates the configuration 
%in the vertical direction while assuming that the instrument points vertically.
%Hence, it corrects for the 25 degree slant angle in the nominal case.

for i=1:12
   dum = fgetl(fid);
end
system17 = str2num(dum(initial_col:end));
Trecv = system17/(256*111.1e3); %in us

% if sver > 1.03              % starting with version 1.05B3,cell size is rounded, not truncated
% 
%    if bin_size == 0.03
%       bin_size = 0.026;
%    elseif bin_size == 0.01
%       bin_size = 0.013;
%    elseif bin_size == 0.04
%       bin_size = 0.039;
%    elseif bin_size == 30
%       %bin_size = 0.026;
%       bin_size = 0.03;
%    elseif bin_size == 40
%       bin_size = 0.039;
%    elseif bin_size == 10
%       bin_size = 0.013;
%    end
% 
% 
% else
%    if bin_size == 0.01
%       bin_size = 0.013;
%    elseif bin_size == 0.02
%       bin_size = 0.026;
%    elseif bin_size == 0.03
%       bin_size = 0.039;
%    elseif bin_size == 20
%       bin_size = 0.02;
%    elseif bin_size == 30
%       bin_size = 0.039;
%    elseif bin_size == 10
%       bin_size = 0.013;
%    end
% end
while ~feof(fid)
   dum = fgetl(fid);
   if length(dum > 0) & strcmp(dum(1:9),'Head freq')
       head_freq = str2num( dum(initial_col:initial_col+4) )*1000;
       break;
   end
end

Tform = NaN*ones(3,3);

while ~feof(fid)
   dum = fgetl(fid);
   if length(dum > 0) & strcmp(dum(1:5),'Trans')
      break;
      indx = find(dum == ' '); 
      indx1 = find(indx < length('Transformation matrix'));
      indx(indx1) = [];
      indx1 = find(diff(indx) > 1);
      Tform(1,1) = str2num(dum(indx(indx1(1))+1:indx(indx1(2))) );
      Tform(1,2) = str2num(dum(indx(indx1(2))+1:indx(end) ));
      Tform(1,3) = str2num(dum(indx(end)+1:end)); 

      for j = 2:3
         dum = fgetl(fid);
         indx = find(dum == ' ');                         
         indx1 = find(indx < length('Transformation matrix'));
         indx(indx1) = [];
         indx1 = find(diff(indx) > 1);
         Tform(j,1) = str2num(dum(indx(indx1(1))+1:indx(indx1(2))) );
         Tform(j,2) = str2num(dum(indx(indx1(2))+1:indx(end) ));
         Tform(j,3) = str2num(dum(indx(end)+1:end)); 
      end

   end

end

fclose(fid);

%NOT SURE WHY THIS IS HARD CODED HERE - SFZ
%for falkor data remove later
%bin_size=0.03;

Head.nmeas = nmeas;
Head.sample_rate= sample_rate;
Head.bin_size = bin_size;
Head.blank_dist = blank_dist;
nblank = ceil(blank_dist/bin_size);
Head.nblank = nblank;
Head.ncells = ncells;
Head.burst_interval = burst_interval;
Head.samples_per_burst = samples_per_burst;
Head.nbeams = nbeams;
Head.pings_per_sample = pings_per_burst;
Head.software_version = sver;
Head.Txmit = Txmit;
Head.Trecv = Trecv;
Head.freq = head_freq;

if ~isnan(Tform(1,1))
   Head.Transformation_matrix = Tform;
end

