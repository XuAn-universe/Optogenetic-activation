function EEG_fdomain_analysis(signal, SampleRate, window_duration, interval)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
Fs = SampleRate;                        % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
cutoff_f = 130;                          % default 130 Hz
% analyze the whole signal
signal = reshape(signal, length(signal), 1);
eeg = signal.*hanning(length(signal));  % signal * window
NFFT0 = 2^nextpow2(length(eeg));        % Next highest power of 2 greater than length(x).
FFTX = fft(eeg, NFFT0);                 % Take FFT, padding with zeros. length(FFTX)==NFFT
NumUniquePts = ceil((NFFT0+1)/2);
FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
MX = MX/length(eeg);                    % Scale the FFT so that it is not a function of the length of x.
MX = MX.^2;                             % this is the power of the signal
f = (0:NumUniquePts-1)*2*Fn/NFFT0;      % frequency coordinate
percentp = MX/sum(MX)*100;              % percentage of power of different frequencies
figure;
hp = plot(f(f <= cutoff_f), percentp(f <= cutoff_f), '-k'); % only plot frequencies less than 130 Hz (you can change this)
yrange = max(percentp(f <= cutoff_f))-min(percentp(f <= cutoff_f)); % the range of the y axis
ylim = [min(percentp(f <= cutoff_f))-0.05*yrange max(percentp(f <= cutoff_f))+0.05*yrange]; % adjust the limits of the y axis
set(gca, 'YLim', ylim, 'XLim', [0 cutoff_f]);
box off;
xlabel('Frequency (Hz)');
ylabel('Percentage of Power (%)');
legend(hp, 'Percentage of Power');
% define different frequency ranges
delta =(f < 4);
theta = (f >= 4 & f < 8);
alpha = (f >= 8 & f < 14);
beta = (f >= 14 & f < 30);
gama = (f >= 30 & f < 130);
% compute the percentage of power in different ranges
percentb(1) = sum(MX(delta))/sum(MX);
percentb(2) = sum(MX(theta))/sum(MX);
percentb(3) = sum(MX(alpha))/sum(MX);
percentb(4) = sum(MX(beta))/sum(MX);
percentb(5) = sum(MX(gama))/sum(MX);
figure
bar(percentb, 'FaceColor', 'black', 'EdgeColor', 'white', 'BarWidth', 0.8);
box off;
set(gca, 'XTick', [1:5], 'XTickLabel', {'Delta', 'Theta', 'Alpha', 'Beta', 'Gama'}, 'TickDir', 'out');
ylabel('Percentage (%)');

% analyze the power spectrum of the signal
wduration = round(window_duration*SampleRate); % transform window duration (in seconds) to data points
dinterval = round(interval/1000*SampleRate); % transform the interval (in miliseconds) between two windows into data points
winfun = hanning(wduration); % just baidu hanning window to see why use it here
NFFT = 2^nextpow2(wduration);          % Next highest power of 2 greater than window duration.
NumUniquePts = ceil((NFFT+1)/2);
f = (0:NumUniquePts-1)*2*Fn/NFFT;
fmap = []; % frequency analysis map
n = 0;
for j = 1:dinterval:length(signal)-wduration+1
    eeg = signal(j:j+wduration-1);          % segment the whole eeg signal to the window duration you defined
    eeg = eeg.*winfun;
    FFTX = fft(eeg, NFFT);                  % Take FFT, padding with zeros. length(FFTX)==NFFT
    FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
    MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
    MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
    MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
    MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
    MX = MX/length(eeg);                    % Scale the FFT so that it is not a function of the length of x.
    MX = MX.^2;                             % analyze the power
    fmap(:, n+1) = MX(f <= cutoff_f);            % keep frequencies below 130 Hz
    n = n+1;
end
x = (1:1:length(signal))/SampleRate; % time of the signal points
xtime = x(1+floor(wduration/2):dinterval:length(signal)-wduration+1+floor(wduration/2))'; % the time of different windows
yfrequency = f(f <= cutoff_f)';
figure;
imshow(flipud(10*log10(fmap)), []); % transform the power to decibel (dB)
colormap jet;
hcb = colorbar;
title(hcb, 'Power (dB)');
axis on;
set(gca, 'XTick', 1:floor(numel(xtime)/6):numel(xtime), 'XTickLabel', num2str(xtime([1:floor(numel(xtime)/6):numel(xtime)]), '%.1f')); % show ~6 ticks for the time axis
set(gca, 'YTick', fliplr(numel(yfrequency)+1-[1:floor(numel(yfrequency)/6):numel(yfrequency)]), 'YTickLabel', num2str(flipud(yfrequency([1:floor(numel(yfrequency)/6):numel(yfrequency)])), '%.1f'));
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Power Spectrum');

