%%
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June. 2015
% xan@cshl.edu
%*---------------------------------------------------------------------*
voltage = 3.3;
Lsensitivity = 0.66;
Lg = (y(:, 1:3)-voltage/2)/Lsensitivity;
Lg = -Lg;

Lxphi = atan2d(Lg(:, 2), Lg(:, 3));
Lytheta = atand(-Lg(:, 1)./sqrt(Lg(:, 2).^2+Lg(:, 3).^2));
figure;
hp = plot(x, Lxphi, '-b');
hold on;
yrange = max(Lxphi)-min(Lxphi);
ylim = [min(Lxphi)-0.05*yrange max(Lxphi)+0.05*yrange];
plot([APre APre], ylim, ':k');
plot([ADuration-APost ADuration-APost], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Angel (degree)');
legend(hp, 'phi');

figure;
hp = plot(x, Lytheta, '-b');
hold on;
yrange = max(Lytheta)-min(Lytheta);
ylim = [min(Lytheta)-0.05*yrange max(Lytheta)+0.05*yrange];
plot([APre APre], ylim, ':k');
plot([ADuration-APost ADuration-APost], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Angel (degree)');
legend(hp, 'theta');

Lxphi0 = Lxphi(APre*SampleRate);
Lytheta0 = Lytheta(APre*SampleRate);
Lgstandard = zeros(size(Lg));
for i = 1:length(x)
    Lgstandard(i, :) = ([cosd(Lytheta0) 0 -sind(Lytheta0); 0 1 0; sind(Lytheta0) 0 cosd(Lytheta0)]^(-1)*[1 0 0; 0 cosd(Lxphi0) sind(Lxphi0); 0 -sind(Lxphi0) cosd(Lxphi0)]^(-1)*[Lg(i, 1); Lg(i, 2); Lg(i, 3)])';
end
figure;
hp = plot(x, Lgstandard(:, 1), '-r', x, Lgstandard(:, 2), '-g', x, Lgstandard(:, 3), '-b');
hold on;
yrange = max(max(Lgstandard)) - min(min(Lgstandard));
ylim = [min(min(Lgstandard))-0.05*yrange max(max(Lgstandard))+0.05*yrange];
plot([APre APre], ylim, ':k');
plot([ADuration-APost ADuration-APost], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Acceleration (g)');
legend(hp, 'x', 'y', 'z');

Lamplitude = sqrt(sum(Lg.^2, 2));
Lthreshold = 5*std(Lamplitude(1:APre*SampleRate));
Lmean = mean(Lamplitude(1:APre*SampleRate));
temp = find(Lamplitude > Lmean+Lthreshold);
temp = temp(temp > APre*SampleRate);
if ~isempty(temp)
    Lstartpointp = temp(1);
else
    Lstartpointp = [];
end
temp = find(Lamplitude < Lmean-Lthreshold);
temp = temp(temp > APre*SampleRate);
if ~isempty(temp)
    Lstartpointn = temp(1);
else
    Lstartpointn = [];
end
if ~isempty(Lstartpointp) && ~isempty(Lstartpointn)
    Lstartpoint = min(Lstartpointp, Lstartpointn);
elseif ~isempty(Lstartpointp) && isempty(Lstartpointn)
    Lstartpoint = Lstartpointp;
elseif isempty(Lstartpointp) && ~isempty(Lstartpointn)
    Lstartpoint = Lstartpointn;
else
    Lstartpoint = [];
end

figure;
hp = plot(x, Lamplitude, '-b');
hold on;
plot([x(1) x(end)], [Lmean+Lthreshold Lmean+Lthreshold], ':k');
plot([x(1) x(end)], [Lmean-Lthreshold Lmean-Lthreshold], ':k');
yrange = max(Lamplitude) - min(Lamplitude);
ylim = [min(Lamplitude)-0.05*yrange max(Lamplitude)+0.05*yrange];
plot([APre APre], ylim, ':k');
plot([ADuration-APost ADuration-APost], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Acceleration (g)');
legend(hp, 'Amplitude');
if ~isempty(Lstartpoint)
    title(['Movement starts at ' num2str((x(Lstartpoint)-APre)*1000) ' ms']);
end

Intperiod = 30;
if ~isempty(Lstartpoint)
    Lgx = Lgstandard(x>=x(Lstartpoint)&x<=x(Lstartpoint)+Intperiod/1000, 1);
    Lgy = Lgstandard(x>=x(Lstartpoint)&x<=x(Lstartpoint)+Intperiod/1000, 2);
    Lgz = Lgstandard(x>=x(Lstartpoint)&x<=x(Lstartpoint)+Intperiod/1000, 3);
    figure;
    hp1 = plot(Lgx, '-r');
    hold on;
    hp2 = plot(Lgy, '-g');
    hp3 = plot(Lgz, '-b');
    yrange = max([Lgx; Lgy; Lgz]) - min([Lgx; Lgy; Lgz]);
    ylim = [min([Lgx; Lgy; Lgz])-0.05*yrange max([Lgx; Lgy; Lgz])+0.05*yrange];
    set(gca, 'YLim', ylim);
    box off;
    ylabel('Acceleration (g)');
    legend([hp1 hp2 hp3], 'x', 'y', 'z');
    Lgx = sum(Lgx);
    Lgy = sum(Lgy);
    Lgz = sum(Lgz);
    ahorizontal = atan2d(Lgy, Lgx);
    aelevation = atand(Lgz/sqrt(Lgx^2+Lgy^2));
    title(['The angle of horizontal is ' num2str(ahorizontal) ' degree and elevation is ' num2str(aelevation) ' degree']);
end

signal = Lamplitude(1:APre*SampleRate);
Fs = SampleRate;                % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
NFFT = 2^nextpow2(length(signal));          % Next highest power of 2 greater than length(x).
FFTX = fft(signal, NFFT);                    % Take FFT, padding with zeros. length(FFTX)==NFFT
NumUniquePts = ceil((NFFT+1)/2);
FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
MX = MX/length(signal);                     % Scale the FFT so that it is not a function of the length of x.
MX = MX.^2;                                        % Power Spectrum
f = (0:NumUniquePts-1)*2*Fn/NFFT;
for i = 1:length(MX)
    if sum(MX(1:i))/sum(MX) >= 0.8
        break;
    end
end
MX = 10*log10(MX);                          % dB
yrange = max(MX(1:i)) - min(MX(1:i));
figure;
plot(f(1:i), MX(1:i), '-k');
xlabel('Frequency (Hz)');
ylabel('Power Spectrum Magnitude (dB)');
set(gca, 'XLim', [-2*Fn/NFFT f(i)+2*Fn/NFFT], 'YLim', [min(MX(1:i))-0.05*yrange max(MX(1:i))+0.05*yrange]);

signal = Lamplitude(APre*SampleRate+1:(ADuration-APost)*SampleRate);
Fs = SampleRate;                % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
NFFT = 2^nextpow2(length(signal));          % Next highest power of 2 greater than length(x).
FFTX = fft(signal, NFFT);                    % Take FFT, padding with zeros. length(FFTX)==NFFT
NumUniquePts = ceil((NFFT+1)/2);
FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
MX = MX/length(signal);                     % Scale the FFT so that it is not a function of the length of x.
MX = MX.^2;                                        % Power Spectrum
f = (0:NumUniquePts-1)*2*Fn/NFFT;
for i = 1:length(MX)
    if sum(MX(1:i))/sum(MX) >= 0.8
        break;
    end
end
MX = 10*log10(MX);                          % dB
yrange = max(MX(1:i)) - min(MX(1:i));
figure;
plot(f(1:i), MX(1:i), '-k');
hold on;
if SP.Frequency == 60
    SP.Frequency = 120;
end
plot([SP.Frequency SP.Frequency], [min(MX(1:i)) max(MX(1:i))], ':k');
xlabel('Frequency (Hz)');
ylabel('Power Spectrum Magnitude (dB)');
set(gca, 'XLim', [-2*Fn/NFFT f(i)+2*Fn/NFFT], 'YLim', [min(MX(1:i))-0.05*yrange max(MX(1:i))+0.05*yrange]);

%%
voltage = 3.3;
Rsensitivity = 0.66;
Rg = (y(:, 4:6)-voltage/2)/Lsensitivity;
Rg = -Rg;

Rxphi = atan2d(Rg(:, 2), Rg(:, 3));
Rytheta = atand(-Rg(:, 1)./sqrt(Rg(:, 2).^2+Rg(:, 3).^2));
figure;
hp = plot(x, Rxphi, '-b');
hold on;
yrange = max(Rxphi)-min(Rxphi);
ylim = [min(Rxphi)-0.05*yrange max(Rxphi)+0.05*yrange];
plot([APre APre], ylim, ':k');
plot([ADuration-APost ADuration-APost], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Angel (degree)');
legend(hp, 'phi');

figure;
hp = plot(x, Rytheta, '-b');
hold on;
yrange = max(Rytheta)-min(Rytheta);
ylim = [min(Rytheta)-0.05*yrange max(Rytheta)+0.05*yrange];
plot([APre APre], ylim, ':k');
plot([ADuration-APost ADuration-APost], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Angel (degree)');
legend(hp, 'theta');

Rxphi0 = Rxphi(APre*SampleRate);
Rytheta0 = Rytheta(APre*SampleRate);
Rgstandard = zeros(size(Rg));
for i = 1:length(x)
    Rgstandard(i, :) = ([cosd(Rytheta0) 0 -sind(Rytheta0); 0 1 0; sind(Rytheta0) 0 cosd(Rytheta0)]^(-1)*[1 0 0; 0 cosd(Rxphi0) sind(Rxphi0); 0 -sind(Rxphi0) cosd(Rxphi0)]^(-1)*[Rg(i, 1); Rg(i, 2); Rg(i, 3)])';
end
figure;
hp = plot(x, Rgstandard(:, 1), '-r', x, Rgstandard(:, 2), '-g', x, Rgstandard(:, 3), '-b');
hold on;
yrange = max(max(Rgstandard)) - min(min(Rgstandard));
ylim = [min(min(Rgstandard))-0.05*yrange max(max(Rgstandard))+0.05*yrange];
plot([APre APre], ylim, ':k');
plot([ADuration-APost ADuration-APost], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Acceleration (g)');
legend(hp, 'x', 'y', 'z');

Ramplitude = sqrt(sum(Rg.^2, 2));
Rthreshold = 5*std(Ramplitude(1:APre*SampleRate));
Rmean = mean(Ramplitude(1:APre*SampleRate));
temp = find(Ramplitude > Rmean+Rthreshold);
temp = temp(temp > APre*SampleRate);
if ~isempty(temp)
    Rstartpointp = temp(1);
else
    Rstartpointp = [];
end
temp = find(Ramplitude < Rmean-Rthreshold);
temp = temp(temp > APre*SampleRate);
if ~isempty(temp)
    Rstartpointn = temp(1);
else
    Rstartpointn = [];
end
if ~isempty(Rstartpointp) && ~isempty(Rstartpointn)
    Rstartpoint = min(Rstartpointp, Rstartpointn);
elseif ~isempty(Rstartpointp) && isempty(Rstartpointn)
    Rstartpoint = Rstartpointp;
elseif isempty(Rstartpointp) && ~isempty(Rstartpointn)
    Rstartpoint = Rstartpointn;
else
    Rstartpoint = [];
end

figure;
hp = plot(x, Ramplitude, '-b');
hold on;
plot([x(1) x(end)], [Rmean+Rthreshold Rmean+Rthreshold], ':k');
plot([x(1) x(end)], [Rmean-Rthreshold Rmean-Rthreshold], ':k');
yrange = max(Ramplitude) - min(Ramplitude);
ylim = [min(Ramplitude)-0.05*yrange max(Ramplitude)+0.05*yrange];
plot([APre APre], ylim, ':k');
plot([ADuration-APost ADuration-APost], ylim, ':k');
set(gca, 'YLim', ylim);
box off;
xlabel('Time (s)');
ylabel('Acceleration (g)');
legend(hp, 'Amplitude');
if ~isempty(Rstartpoint)
    title(['Movement starts at ' num2str((x(Rstartpoint)-APre)*1000) ' ms']);
end

Intperiod = 30;
if ~isempty(Rstartpoint)
    Rgx = Rgstandard(x>=x(Rstartpoint)&x<=x(Rstartpoint)+Intperiod/1000, 1);
    Rgy = Rgstandard(x>=x(Rstartpoint)&x<=x(Rstartpoint)+Intperiod/1000, 2);
    Rgz = Rgstandard(x>=x(Rstartpoint)&x<=x(Rstartpoint)+Intperiod/1000, 3);
    figure;
    hp1 = plot(Rgx, '-r');
    hold on;
    hp2 = plot(Rgy, '-g');
    hp3 = plot(Rgz, '-b');
    yrange = max([Rgx; Rgy; Rgz]) - min([Rgx; Rgy; Rgz]);
    ylim = [min([Rgx; Rgy; Rgz])-0.05*yrange max([Rgx; Rgy; Rgz])+0.05*yrange];
    set(gca, 'YLim', ylim);
    box off;
    ylabel('Acceleration (g)');
    legend([hp1 hp2 hp3], 'x', 'y', 'z');
    Rgx = sum(Rgx);
    Rgy = sum(Rgy);
    Rgz = sum(Rgz);
    ahorizontal = atan2d(Rgy, Rgx);
    aelevation = atand(Rgz/sqrt(Rgx^2+Rgy^2));
    title(['The angle of horizontal is ' num2str(ahorizontal) ' degree and elevation is ' num2str(aelevation) ' degree']);
end

signal = Ramplitude(1:APre*SampleRate);
Fs = SampleRate;                % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
NFFT = 2^nextpow2(length(signal));          % Next highest power of 2 greater than length(x).
FFTX = fft(signal, NFFT);                    % Take FFT, padding with zeros. length(FFTX)==NFFT
NumUniquePts = ceil((NFFT+1)/2);
FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
MX = MX/length(signal);                     % Scale the FFT so that it is not a function of the length of x.
MX = MX.^2;                                        % Power Spectrum
f = (0:NumUniquePts-1)*2*Fn/NFFT;
for i = 1:length(MX)
    if sum(MX(1:i))/sum(MX) >= 0.8
        break;
    end
end
MX = 10*log10(MX);                          % dB
yrange = max(MX(1:i)) - min(MX(1:i));
figure;
plot(f(1:i), MX(1:i), '-k');
xlabel('Frequency (Hz)');
ylabel('Power Spectrum Magnitude (dB)');
set(gca, 'XLim', [-2*Fn/NFFT f(i)+2*Fn/NFFT], 'YLim', [min(MX(1:i))-0.05*yrange max(MX(1:i))+0.05*yrange]);

signal = Ramplitude(APre*SampleRate+1:(ADuration-APost)*SampleRate);
Fs = SampleRate;                % sampling frequency
Fn = Fs/2;                              % Nyquist frequency
NFFT = 2^nextpow2(length(signal));          % Next highest power of 2 greater than length(x).
FFTX = fft(signal, NFFT);                    % Take FFT, padding with zeros. length(FFTX)==NFFT
NumUniquePts = ceil((NFFT+1)/2);
FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
MX = MX/length(signal);                     % Scale the FFT so that it is not a function of the length of x.
MX = MX.^2;                                        % Power Spectrum
f = (0:NumUniquePts-1)*2*Fn/NFFT;
for i = 1:length(MX)
    if sum(MX(1:i))/sum(MX) >= 0.8
        break;
    end
end
MX = 10*log10(MX);                          % dB
yrange = max(MX(1:i)) - min(MX(1:i));
figure;
plot(f(1:i), MX(1:i), '-k');
hold on;
if SP.Frequency == 60
    SP.Frequency = 120;
end
plot([SP.Frequency SP.Frequency], [min(MX(1:i)) max(MX(1:i))], ':k');
xlabel('Frequency (Hz)');
ylabel('Power Spectrum Magnitude (dB)');
set(gca, 'XLim', [-2*Fn/NFFT f(i)+2*Fn/NFFT], 'YLim', [min(MX(1:i))-0.05*yrange max(MX(1:i))+0.05*yrange]);