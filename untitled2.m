clear;
clc

T = 41.6*1e-6;
N = 16;
alpha_max = 1;

fc = 4e9;
c = 3e8;
lambda = c/fc;

delta = 1/T;
Doppler_freq_max = alpha_max*delta/N;

Doppler_speed_max = Doppler_freq_max * lambda;

fprintf("载波间隔周期：%.2f us.\n", T*1e6)
fprintf("载波频率：%.2f Hz.\n", delta)
fprintf("载波波长：%f.\n", lambda)
fprintf("最大多普勒频率：%f Hz.\n", Doppler_freq_max)
fprintf("多普勒速度 %f m/s\n", Doppler_speed_max);
fprintf("多普勒速度 %f km/h\n\n", Doppler_speed_max*3.6);

NT = N*T;   % 不一定对……
N = 256;
T = NT/N;

alpha_max = 2;

fc = 4e9;
c = 3e8;
lambda = c/fc;

delta = 1/T;
Doppler_freq_max = alpha_max*delta/N;

Doppler_speed_max = Doppler_freq_max * lambda;

fprintf("载波间隔周期：%.2f us.\n", T*1e6)
fprintf("载波频率：%.2f Hz.\n", delta)
fprintf("载波波长：%f.\n", lambda)
fprintf("最大多普勒频率：%f Hz.\n", Doppler_freq_max)
fprintf("多普勒速度 %f m/s\n", Doppler_speed_max);
fprintf("多普勒速度 %f km/h\n", Doppler_speed_max*3.6);



% 
% speed2dop(405/3.6, lambda)
% 
% dop2speed



clear; clc

c = 3e8; fc = 4e9; lambda = c/fc;

% Case 1
Ttap = 41.6e-6; fs = 1/Ttap; N = 16; alpha_max = 1;
fD_max = alpha_max * fs / N;
vmax = (c/fc) * fD_max;
fprintf("Case1: N=%d, αmax=%d, fs=%.1f Hz → v=%.1f km/h\n", N, alpha_max, fs, vmax*3.6);

% Case 2
Ttap = 3.906e-6; fs = 1/Ttap; N = 256; alpha_max = 2;
fD_max = alpha_max * fs / N;
vmax = (c/fc) * fD_max;
fprintf("Case2: N=%d, αmax=%d, fs=%.1f Hz → v=%.1f km/h\n", N, alpha_max, fs, vmax*3.6);




taps = 3;
num_trials = 1e5;

% Case A: normalized total power
hA = 1/sqrt(2*taps) * (randn(num_trials,taps) + 1i*randn(num_trials,taps));
pA = mean(sum(abs(hA).^2,2));

% Case B: each tap unit variance
hB = 1/sqrt(2) * (randn(num_trials,taps) + 1i*randn(num_trials,taps));
pB = mean(sum(abs(hB).^2,2));

fprintf('Average total power (A): %.3f\n', pA);
fprintf('Average total power (B): %.3f\n', pB);
