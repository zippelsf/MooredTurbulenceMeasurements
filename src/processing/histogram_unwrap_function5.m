
function [unwrap_data] = histogram_unwrap_function5(V1_fdata,Vmax,dim)
%
% dim-->dimension to unwrap along

%%For testing new algorithm
%V1_fdata=V1_tdata;
%profile=282;

if dim==2
[numt,SetL]=size(V1_fdata);
elseif dim==1
V1_fdata=V1_fdata';
[numt,SetL]=size(V1_fdata);
end

SetL = length(V1_fdata(1,:));
twopi=2*Vmax;

unwrap_data = V1_fdata;

for profile = 1:numt
clear AA ii nn n mm m bins bins1 Vrep zeros histnums win1 win2 TT binx

if isnan(V1_fdata(profile,:))
   unwrap_data(profile,:) = NaN;
else
        
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

%%%%%%%%%%%%%%%%%%%
% Leah's code discards profiles with 3 or more NaNs.
% This recreates that, though 3 is a bit strict:
% (To exacly recreate, use ">=3")
if length(nan_ind)>=20
  unwrap_data(profile,:) = NaN;
end

end
unwrap_data(profile,:) = unwrap_data(profile,:) -nanmean(unwrap_data(profile,:));
end



