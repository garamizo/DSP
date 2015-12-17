%% DSP - Laboratory #6
% Guilherme Aramizo Ribeiro
%% Definitions
% TODO replace with experimental data
rpm_nominal_ = 2490; % motor operating speed [rpm]
rpm_min_ = 225;
acc_sens_ = 0.03417; % V/(m/s^2)

%% Model of the (t,w) -> acc function

deriv = @(y,Fs) savitzkyGolayFilt(-y,10,1,51) * 1/(1/Fs)^1;
integ = @(y,Fs) cumsum(y)/Fs;
orderFunc = @(w,mr,O,Fs) mr*O*(-w.^2.*sin(O*integ(w,Fs))*O + deriv(w,Fs).*cos(O*integ(w,Fs))); % force out

sys = tf([1 0 0], [5 1 3e3]); % mass spring damper system
motor_func_ = @(w,Fs) lsim(sys, orderFunc(w, 3e-6, 1, Fs) + ...
    orderFunc(w, 2e-6, 2, Fs) + orderFunc(w, 1e-6, 5.3, Fs) + ...
    orderFunc(w, 0.5e-6, 8, Fs) + orderFunc(w, 1e-6, 10, Fs), (0:length(w)-1)/Fs); % acc out

% motor_func_ = @(w,Fs) lsim(sys, orderFunc(w, 1e-2, 1, Fs), (0:length(w)-1)/Fs); % acc out

%{
Fs = 13.1072e6/(256*5);
time = 0 : (1/Fs) : 10;
w = linspace(rpm_min_, rpm_nominal_, length(time))*2*pi/60;

plot(time, motor_func_(time, w))
%}
%% Getting calibration factor
% Calibrating a signal
%     Input frequency range: [Fmin Fmax] = [149.2 149.2] Hz
%     Input voltage range: [Vmin Vmax] = ?
% 
%     Fs = 13.1072e6/(256*31) > 300 Hz = 2*Fmax
%     winsize = 476, minimize leakage
%     noavg = 32

% acquisition parameters
Fs = 13.1072e6/(256*31);
winsize = 476;
noavg = 32;
range = 5;
bits = 24;

% TODO replace with experimental data
% Simulate signal
% time = linspace(0, (winsize*(noavg+2))/Fs, winsize*(noavg+2)).';
% y = myDSP.discretize( ...
%     (9.81*sqrt(2))*acc_sens_*sin(2*pi*149.2*time), bits, range);
data = dlmread('../report_lab06/data/calib.csv', '\t');
time = data(:,1);
y = data(:,2);

%{
plot(time*1e3, y, 'o-')
xlim([0 30])
xlabel('time [ms]'); ylabel('acc voltage [V]')
%}

% Visualize auto-power
yreshaped = myDSP.reshape(y, winsize, 0);

win = window(@flattopwin, winsize);
gain = [1; 2*ones(winsize-1, 1)] / (mean(win)*winsize);

GY = bsxfun(@times, fft(bsxfun(@times, yreshaped, win)), gain);
Gyy = squeeze(mean(conj(GY) .* GY, 2));

f = linspace(0, Fs, winsize);

figure
subplot(211); semilogy(f, sqrt(Gyy)); grid on; grid minor
xlim([0 Fs/2]); xlabel('frequency [Hz]')
subplot(212); plot(time*1e3, y, 'o-'); xlabel('time [s]'); xlim([0 30])

calib_factor = (9.81*sqrt(2))/0.03417; % Calibration factor [(m/s^2)/V]

%% Set #1

% Steady state 0.5*rpm speed
%     Input frequency range: 10th order of the nominal speed
%     Input voltage range: [Vmin Vmax] = ?
% 
%     Fs = 51200 = 13.1072e6/(256*1) > 2*Fmax = 2*(10*(rpm_nominal_/60))
%     winsize = 204800 = Fs/0.25, df = Fs/winsize = 0.25
%     noavg = 50

% acquisition parameters
Fs = 13.1072e6/(256*31);% > 2*Fmax = 2*(10*(rpm_nominal_/60))
winsize = 6606; % = Fs/0.25, df = Fs/winsize = 0.25
noavg = 50;
range = 5;
bits = 24;

% TODO replace with experimental data
time = linspace(0, (noavg+2)*winsize/Fs, (noavg+2)*winsize).';
w = ones(size(time)) *(rpm_nominal_/2)*2*pi/60;
tach = square(cumsum(w)/Fs) + randn(size(time))/10;
acc = calib_factor * myDSP.discretize( ...
    acc_sens_ * motor_func_(w, Fs), bits, range);

% data = dlmread('../report_lab06/data/ss.csv', '\t', 1);
% time = data(:,1) - data(1,1);
% tach = data(:,4);
% acc = data(:,2) * calib_factor;

%{
figure
ax2 = subplot(211); stairs(time, tach); xlim([0 180e-3])
ax3 = subplot(212); stairs(time, acc); xlim([0 180e-3])
linkaxes([ax2,ax3],'x')
%}

yreshaped = myDSP.reshape(acc, winsize, 0);
yreshaped(:,[1 end]) = []; % remove head and tail

win = window(@hann, winsize);

gain = [1; 2*ones(winsize-1, 1)] / (mean(win)*winsize);

GY = bsxfun(@times, fft(bsxfun(@times, yreshaped, win)), gain);
Gyy = squeeze(mean(conj(GY) .* GY, 2));

f = linspace(0, Fs, winsize);

figure
semilogy(f, sqrt(Gyy)); grid on; grid minor
xlim([0 Fs/2]); xlabel('frequency [Hz]'); ylabel('base acceleration [m/s^2]')

%% 1)
avg_rpm = mean(myDSP.speed_from_tach(time, tach)); % Average RPM

%%
idx = time > 5 & time < 10;
myDSP.TimeTVDFTOrderTracker(acc(idx), Fs, w(idx)*60/(2*pi));

%% Set #2, #3, #4

% Transient: 0->45s, 0->20s, 0->10s, ranging full rpm bounds
%     Input frequency range: 10th order of nominal speed
%     Input voltage range: [Vmin Vmax] = ?
% 
%     Fs = 10240 = 13.1072e6/(256*5) > 2*Fmax = 2*(10*(rpm_nominal_/60))
%     winsize = 40960 = Fs/0.25, df = Fs/winsize = 0.25
%     noavg = 50

% acquisition parameters
Fs = 13.1072e6/(256*31);
len = 45*Fs;
% winsize = 40960;
% noavg = 50;
range = 5;
bits = 24;

time = linspace(0, round(len*1.1)/Fs, round(len*1.1)).';

w = ones(size(time, 1), 3) * rpm_nominal_ * 2*pi/60;
w(time < 45, 1) = (rpm_min_ + time(time < 45)*(rpm_nominal_-rpm_min_)/45)*2*pi/60;
w(time < 20, 2) = (rpm_min_ + time(time < 20)*(rpm_nominal_-rpm_min_)/20)*2*pi/60;
w(time < 10, 3) = (rpm_min_ + time(time < 10)*(rpm_nominal_-rpm_min_)/10)*2*pi/60;

w(:,1) = w(:,1) - 30*sin(pi*time/time(end));

th = cumsum(w, 1)/Fs;
tach = square(th) + randn(size(th))/10;

acc = calib_factor * myDSP.discretize(acc_sens_ * [motor_func_(w(:,1), Fs), ...
    motor_func_(w(:,2), Fs), motor_func_(w(:,3), Fs)], bits, range);

%{
data1 = dlmread('../report_lab06/data/slow3.csv', '\t', 1);
data2 = dlmread('../report_lab06/data/mid2.csv', '\t', 1);
data3 = dlmread('../report_lab06/data/fast2.csv', '\t', 1);

len = max([size(data1,1), size(data2,1), size(data3,1)]);
time = (1:len).'/Fs;

tach = zeros(len, 3);
acc = zeros(len, 3);

tach(1:size(data1,1),1) = data1(:,4);
tach(1:size(data2,1),2) = data2(:,4);
tach(1:size(data3,1),3) = data3(:,4);

acc(1:size(data1,1),1) = data1(:,2) * calib_factor;
acc(1:size(data2,1),2) = data2(:,2) * calib_factor;
acc(1:size(data3,1),3) = data3(:,2) * calib_factor;
%}

%% 2) Tachometer
rpm = [myDSP.speed_from_tach(time, tach(:,1)), myDSP.speed_from_tach( ...
    time, tach(:,2)), myDSP.speed_from_tach(time, tach(:,3))];

figure
plot(time, rpm); grid on; grid minor
xlabel('frequency [Hz]'); ylabel('Shaft speed [RPM]'); ylim([rpm_min_*0.9 rpm_nominal_*1.1])

%% 3) Color map

t0 = [0 0 0];
tend = [45 20 10];
for k = 1 : 3
    idx = time > t0(k) & time < tend(k);
    len = sum(idx); % # samples of test
    winsize = floor(len/100); % 100 evenly spaced RPM bins

    accblock = myDSP.reshape(acc(idx, k), winsize, 0);
    rpmblock = myDSP.reshape(rpm(idx,k), winsize, 0);

    win = window(@hann, winsize);
    gain = [1; 2*ones(winsize-1, 1)] / (rms(win)*winsize);

    GY = bsxfun(@times, fft(bsxfun(@times, accblock, win)), gain);
    Gyy = conj(GY) .* GY;

    f = linspace(0, Fs, winsize);
    rpmOfBlock = mean(rpmblock, 1);
    rpmSpaced = linspace(min(rpmOfBlock), max(rpmOfBlock), length(rpmOfBlock));
    [F, RPM] = meshgrid(f, rpmSpaced);
    accs = sqrt(Gyy).';
    ACC = interp2(f,rpmOfBlock,accs,F,RPM);

    figure;
    surf(F, RPM, ACC, 'EdgeColor','None');
    view(2); xlim([0 10*rpm_nominal_/60]); ylim([rpm_min_, rpm_nominal_])
    colorbar; xlabel('Frequency [Hz]'); ylabel('shaft speed [RPM]')
end

%% 4) Order tracking
%% FFT, constant order bandwidth

myDSP.colorMap(o, r, g) 

rpm = w * 60/(2*pi);
idx = time > 5 & time < 45;

for k = 1 : 3
    [o, r, g] = myDSP.OrderResamplingOrderTracker(acc(idx,2), Fs, rpm(idx,2));
    [o, r, g] = myDSP.OrderTVDFTOrderTracker(acc(idx,1), Fs, rpm(idx,1));
    [o, r, g] = myDSP.TimeTVDFTOrderTracker(acc(idx,1), Fs, rpm(idx,1));


