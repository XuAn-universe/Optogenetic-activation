function [f, MX] = FFT_trajectory(trajectory, Fn)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
winfunc = hanning(numel(trajectory));
NFFT = 2^nextpow2(length(winfunc));         % Next highest power of 2 greater than length(x).
NumUniquePts = ceil((NFFT+1)/2);
f = (0:NumUniquePts-1)*2*Fn/NFFT;
signal = (trajectory-mean(trajectory)).*winfunc;
FFTX = fft(signal, NFFT);               % Take FFT, padding with zeros. length(FFTX)==NFFT
FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
MX = abs(FFTX).^2;                       
MX(2:end-1) = 2*MX(2:end-1);
MX = MX/sum(winfunc.^2);
