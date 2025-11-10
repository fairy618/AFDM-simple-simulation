%% FOR Fig. 10. BER and spectral efficiency performance of OFDM, OCDM, OTFS and AFDM using MMSE detection.

clear; clc;
% close all;
% rng(7)
tic
% System parameters
M_mod = 4;
N = 256;
T_tap = 3.9e-6;
N_frame = 1e3; 
SNR_dB = 0:2:20;  

fc = 4e9;
c = 3e8;
lambda = c/fc;

fs = 1/T_tap;

T_sym = N * T_tap;
delta_f = 1/T_sym;

eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));   % average power per symbol

SNR = 10.^(SNR_dB/10);
k = log2(M_mod);
sigma_2 = (abs(eng_sqrt)^2)./(SNR * k);

% Generate synthetic delay-Doppler channel
k_max = 2;  % maximum normalized Doppler index
taps  = 3;  % number of paths
l_max = 2;  % maximum normalized delay index

chan_coef = 1/sqrt(2*taps).*(randn(1,taps)+1i.*randn(1,taps));  % 总能量固定

delay_taps = randi(l_max, [1,taps]) - 1;

fD_max = k_max / T_sym;  % 最大物理多普勒频移
theta = (rand(1,taps)*2 - 1) * pi;   % uniform in [-pi, pi]
Doppler_freq = fD_max * cos(theta);         % Hz
Doppler_taps = Doppler_freq * T_sym;

trellis = poly2trellis(7,[171 133]);

% -------- Summary printout --------
v_max = (c/fc) * fD_max;    % 最大多普勒对应速度 (m/s)
v_i   = (c/fc) * Doppler_freq;   % 每条路径的速度 (m/s)

fprintf('\n=========== System Summary ===========\n');
fprintf('Carrier frequency (fc):       %.3f GHz\n', fc/1e9);
fprintf('Wavelength (lambda):          %.4f m\n', lambda);
fprintf('Delay-tap spacing (T_tap):    %.3f µs\n', T_tap*1e6);
fprintf('Sampling rate (fs):           %.3f kHz\n', fs/1e3);
fprintf('Symbol duration (T_sym):      %.3f µs\n', T_sym*1e6);
fprintf('Subcarrier spacing (Δf):      %.3f Hz\n', delta_f);
fprintf('Grid size (N):                %d\n', N);
fprintf('Max. normalized Doppler (k_max): %d\n', k_max);
fprintf('Max. physical Doppler (fD_max): %.3f Hz\n', fD_max);
fprintf('Max. Doppler speed (v_max):   %.3f m/s  (%.3f km/h)\n', v_max, v_max*3.6);
fprintf('Number of paths (taps):       %d\n', taps);
fprintf('Max. normalized delay index (l_max): %d\n', l_max);

fprintf('\n---- Per-path Doppler info ----\n');
for i = 1:taps
    fprintf('Path %d: delay=%d, Doppler=%.3f Hz, speed=%.3f m/s (%.2f km/h), Doppler_tap=%.3f\n', ...
        i, delay_taps(i), Doppler_freq(i), v_i(i), v_i(i)*3.6, Doppler_taps(i));
end
fprintf('================================\n\n');


% AFDM parameters %%
max_Doppler = max(Doppler_taps);
max_delay = max(delay_taps);

CPP_len = max_delay;    % CPP_len >= l_max-1
N_data = N-CPP_len;     % length of data symbols

CP_len = ceil(max(delay_taps)) + 2;
N_data_ofdm = N - CP_len;

k_v = 1;    % guard interval to combat fractional Doppler shifts, see equation (38) in [R1]
if (2*(max_Doppler+k_v)*(max_delay+1)+max_delay)>N_data
    error('subcarrier orthogonality is not satisfied');
end
c1 = (2*(max_Doppler+k_v)+1)/(2*N_data);    % equation (48) in [R1]
c2 = 1/(N_data^2);

fprintf("max_Doppler = %.2f. max_delay = %d.\n\n", max_Doppler, max_delay);

% Generate channel matrix
L_set = unique(delay_taps);
q = 0:N-1; % 所有频率索引
phase = exp(-1i*2*pi*(Doppler_freq(:) * q));   % taps × N
weighted = chan_coef(:) .* phase;  % taps × N
gs=zeros(max_delay+1,N);
for i = 1:taps
    gs(delay_taps(i)+1, :) = gs(delay_taps(i)+1, :) + weighted(i, :);
end
gsConst = parallel.pool.Constant(gs);

% channel matrix form
H = Gen_channel_mtx(N, taps, chan_coef, delay_taps, Doppler_freq, c1);  % equation (24) in [R1]
% Observe the structure of H
% imagesc(abs(H))

%% Start Loop
ber_AFDM  = zeros(size(SNR_dB));

for iesn0 = 1:length(SNR_dB)

    sigma2 = sigma_2(iesn0);
    err_count_AFDM = zeros(N_frame,1);

    parfor iframe = 1:N_frame
        % Tx data generation %%
        x = randi([0, M_mod-1], N_data, 1);
        x_qam = qammod(x, M_mod, 'gray', 'UnitAveragePower', true);
        w = sqrt(sigma2/2) * (randn(N, 1) + 1i*randn(N, 1));

        % AFDM chain
        s_afdm = AFDM_mod(x_qam, c1, c2);
        cpp_afdm = s_afdm(N_data-CPP_len:N_data-1).*exp(-1i*2*pi*c1*(N^2+2*N*(-CPP_len:-1).'));
        s_cpp_afdm = [cpp_afdm; s_afdm];

        gs_local = gsConst.Value;
        r_afdm = zeros(N,1);
        for l = (L_set+1)
            r_afdm(l:N) = r_afdm(l:N) + gs_local(l, l:N).' .* s_cpp_afdm(1:N-l+1);
        end
        r_afdm = r_afdm + w;

        x_est_afdm = H'/(H*H'+sigma2*eye(N))*r_afdm;
        x_est_no_cpp_afdm = x_est_afdm(CPP_len+1:end);
        y_afdm = AFDM_demod(x_est_no_cpp_afdm, c1, c2);

        x_est_bit_afdm = qamdemod(y_afdm, M_mod, 'gray', 'UnitAveragePower', true);

        err_count_AFDM(iframe) = sum(x ~= x_est_bit_afdm);

    end
    ber_AFDM(iesn0) = sum(err_count_AFDM)/(N_data * N_frame);
    fprintf('SNR=%2d dB done: AFDM=%.3e\n', SNR_dB(iesn0), ber_AFDM(iesn0));
end

% Plot bit error rate
figure;
semilogy(SNR_dB, ber_AFDM, '-o', 'LineWidth', 1.2);
grid on;
xlabel('SNR (dB)');
ylabel('Symbol error rate (per symbol)');
title(sprintf("N=%d P=%d K=%d L=%d", N, taps, k_max, l_max));

toc

