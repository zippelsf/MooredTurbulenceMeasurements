%[yfilt]=synthestic_timeseries(N,beta);
%
%Simulate a time series with a given spectral power law
%
%Inputs:
%  N, number of points
%  Beta, power-law dependence, such that spectrum is proportional to freq^-Beta
%
%Output:
%  yfilt, a unit-variance timerseries with the specified power-law dependence
%  
% Note that "unit variance" is imposed by brute force and is exact, not statistical.
% That is, the computed variance will be exactly 1.
%
% Started 12/28/2015
% Tom Farrar, jfarrar@whoi.edu
%  Maybe in a later version I should allow freq and spectrum to be input, along with N and dt
%



function [yfilt]=synthetic_timeseries(N,beta);

%%%%%%%%%%%%%
% Test values:
%N=1000;
%beta=-2;
%%%%%%%%%%%%%

dt=1;
y=randn(N,1);
y=y(:);

%Generate frequency axis
if mod(N,2)==0
    k=-N/2:N/2-1; % N even
else
    k=-(N-1)/2:(N-1)/2; % N odd
end
T=N.*dt;
freq=k/T;  %the frequency axis
freq=freq(:);

%takes the fft of the signal, and adjusts the amplitude accordingly
Y=fft(y); % 
Y=fftshift(Y); %shifts the fft data so that it is centered and aligns with freq

Yfilt=abs(freq).^(beta./2).*Y;
ff=find(freq==0);Yfilt(ff)=0;
Yfilt2=ifftshift(Yfilt);

yfilt=ifft(Yfilt2);

%Lines below confirm that fftshift, fft, and inverses all work as expected and give original noise timeseries:
%Ytest=ifftshift(Y);
%ytest=ifft(Ytest);
%max_ytest_error=max(y-ytest)

%Normalize to unit variance:
vv=var(yfilt);
yfilt=yfilt./sqrt(vv);
