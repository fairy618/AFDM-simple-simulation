clear; clc;
tic

M = 4;              % 调制
N = 256;            % 子载波数量
EbN0_dB = 0:2:20;   % 信噪比

T_tap = 3.9e-6;  

N_frame = 1e2;      % 仿真次数

fc = 4e9;
c = 3e8;
lambda = c/fc;

fs = 1/T_tap;
T_sym = N * T_tap;
delta_f = 1/T_sym;

eng_sqrt = (M==2)+(M~=2)*sqrt((M-1)/6*(2^2));   % average power per symbol

k = log2(M);
EbN0 = 10.^(EbN0_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./(EbN0 * k);

k_max = 2;  % maximum normalized Doppler index
taps  = 3;  % number of paths
l_max = 2;  % maximum normalized delay index

chan_coef = 1/sqrt(2*taps).*(randn(1,taps)+1i.*randn(1,taps));  % 总能量固定

delay_taps = randi(l_max, [1,taps]) - 1;

fD_max = k_max / T_sym;  % 最大物理多普勒频移
theta = (rand(1,taps)*2 - 1) * pi;   % uniform in [-pi, pi]
Doppler_freq = fD_max * cos(theta);         % Hz
Doppler_taps = Doppler_freq * T_sym;

max_Doppler = max(Doppler_taps);
max_delay = max(delay_taps);

CPP_len = max_delay;    % CPP_len >= l_max-1
N_data = N-CPP_len;     % length of data symbols

Ncp = round(0.07 * N);  % CP 长度 (取整)
ofdmSym = 1;            % 每次发送多少 OFDM 符号
N_data_ofdm = N - Ncp;

k_v = 1;    % guard interval to combat fractional Doppler shifts, see equation (38) in [R1]
if (2*(max_Doppler+k_v)*(max_delay+1)+max_delay)>N_data
    error('subcarrier orthogonality is not satisfied');
end
c1 = (2*(max_Doppler+k_v)+1)/(2*N_data);    % equation (48) in [R1]
c2 = 1/(N_data^2);

fprintf("max_Doppler = %.2f. max_delay = %d.\n\n", max_Doppler, max_delay);

L_set = unique(delay_taps);
q = 0:N-1; % 所有频率索引
phase = exp(-1i*2*pi*(Doppler_freq(:) * q));   % taps × N
weighted = chan_coef(:) .* phase;  % taps × N
gs=zeros(max_delay+1,N);
for i = 1:taps
    gs(delay_taps(i)+1, :) = gs(delay_taps(i)+1, :) + weighted(i, :);
end
gsConst = parallel.pool.Constant(gs);

H = Gen_channel_mtx(N, taps, chan_coef, delay_taps, Doppler_freq, c1);  % equation (24) in [R1]

%% Start Loop
ber_AFDM  = zeros(size(EbN0_dB));
ber_OFDM  = zeros(size(EbN0_dB));

for iesn0 = 1:length(EbN0_dB)

    sigma2 = sigma_2(iesn0);

    err_count_AFDM = zeros(N_frame,1);
    err_count_OFDM = zeros(N_frame,1);

    parfor iframe = 1:N_frame
        % Tx data generation %%
        x = randi([0, M-1], N, 1);
        x_qam = qammod(x, M, 'gray', 'UnitAveragePower', true);
        
        % AFDM chain
        s_afdm = AFDM_mod(x_qam(1:N_data), c1, c2);
        cpp_afdm = s_afdm(N_data-CPP_len:N_data-1).*exp(-1i*2*pi*c1*(N^2+2*N*(-CPP_len:-1).'));
        s_cpp_afdm = [cpp_afdm; s_afdm];

        gs_local = gsConst.Value;
        r_afdm = zeros(N,1);
        for l = (L_set+1)
            r_afdm(l:N) = r_afdm(l:N) + gs_local(l, l:N).' .* s_cpp_afdm(1:N-l+1);
        end
        w_afdm = sqrt(sigma2 / (2)) * (randn(N, 1) + 1i*randn(N, 1));
        r_afdm = r_afdm + w_afdm;
        x_est_afdm = H'/(H*H'+sigma2*eye(N))*r_afdm;
        
        x_est_no_cpp_afdm = x_est_afdm(CPP_len+1:end);
        y_afdm = AFDM_demod(x_est_no_cpp_afdm, c1, c2);

        x_est_bit_afdm = qamdemod(y_afdm, M, 'gray', 'UnitAveragePower', true);

        err_count_AFDM(iframe) = sum(x(1:N_data) ~= x_est_bit_afdm);

        % OFDM chain
        s_ofdm = ifft(x_qam(1:N_data_ofdm), N_data_ofdm, 1);
        cp_ofdm = s_ofdm(end-Ncp+1:end);
        s_cp_ofdm = [cp_ofdm; s_ofdm];
        % 待添加信道
        w_ofdm = sqrt(sigma2 / (2*N_data_ofdm)) * (randn(N, 1) + 1i*randn(N, 1));
        r_ofdm = s_cp_ofdm + w_ofdm;
        x_est_ofdm = r_ofdm;
        x_est_no_cpp_ofdm = x_est_ofdm(Ncp+1:end);
        y_ofdm = fft(x_est_no_cpp_ofdm, N_data_ofdm, 1);

        x_est_bit_ofdm = qamdemod(y_ofdm, M, 'gray', 'UnitAveragePower', true);

        err_count_OFDM(iframe) = sum(x(1:N_data_ofdm) ~= x_est_bit_ofdm);
    end
    ber_AFDM(iesn0) = sum(err_count_AFDM)/(N_data * N_frame);
    ber_OFDM(iesn0) = sum(err_count_OFDM)/(N_data_ofdm * N_frame);

    fprintf('Eb/N0=%2d dB done: AFDM=%.3e\tOFDM=%.3e\n', EbN0_dB(iesn0), ber_AFDM(iesn0), ber_OFDM(iesn0));
end

% Plot bit error rate
figure;
semilogy(EbN0_dB, ber_AFDM, '-o', 'LineWidth', 1.2); hold on
semilogy(EbN0_dB, ber_OFDM, '-o', 'LineWidth', 1.2); hold off

grid on;
xlabel('SNR (dB)');
ylabel('Symbol error rate (per symbol)');
legend("AFDM", "OFDM")
title(sprintf("N=%d P=%d K=%d L=%d", N, taps, k_max, l_max));

toc

