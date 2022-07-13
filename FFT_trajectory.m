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
MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
MX = MX/length(signal);                 % Scale the FFT so that it is not a function of the length of x.
MX = MX.^2;