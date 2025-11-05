%% FOR Fig. 10. BER and spectral efficiency performance of OFDM, OCDM, OTFS and AFDM using MMSE detection.

clear;
clc;

tic
% System parameters
M_mod = 4;      % size of QAM constellation
N = 256;        % number of symbols(subcarriers)
B = 10e6;
c = 3e8;

N_frame = 100;    % number of simulation frames

if floor(log2(M_mod)) ~= log2(M_mod)
    error('M_mod must be a power of 2 for bit mapping.');
end

car_fre = 4e9;  % carrier frequency
delta_f = 15e3;  % symbol spacing    符号间距
T = 1/delta_f;  % symbol duration   符号持续时间

k = log2(M_mod);
Rc = 1/2;
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));   % average power per symbol
SNR_dB = -10:2:20;        

SNR = 10.^(SNR_dB/10);  % 符号能量/噪声功率
sigma_2 = (abs(eng_sqrt)^2)./(SNR * k);   % noise power
sigma_2_code = (abs(eng_sqrt)^2) ./ (SNR * k * Rc);

% Generate synthetic delay-Doppler channel %% 生成合成延迟-多普勒信道
k_max = 30;    %  maximum Doppler shift is αmax = 2
taps  = 3;    % a 3-path channel
l_max = 2;    % maximum delay spread

chan_coef = 1/sqrt(2).*(randn(1,taps)+1i.*randn(1,taps));   % follows Rayleigh distribution
fprintf("P = %d. Channel Power = %.2f.\n", taps, sum(abs(chan_coef).^2))

delay_taps = randi(l_max, [1,taps]) - 1;
fprintf("delay_taps:"); disp(delay_taps);

fD_max = k_max / (N*T);  % 最大物理多普勒频移
u = rand(1, taps);
Doppler_freq = fD_max * sin(pi * (u - 0.5));   % 服从近似Jakes分布
Doppler_taps = Doppler_freq * N*T;
% fprintf("Doppler_freq : %.2fkHz.\n", Doppler_freq);
Doppler_speed_ms = Doppler_freq * c / car_fre;
Doppler_speed_kmh = Doppler_speed_ms * 3.6;
fprintf("Doppler_speed_kmh：%.2f km/h.\n", Doppler_speed_kmh)

% AFDM parameters
max_Doppler = max(Doppler_taps);
max_delay = max(delay_taps);

CPP_len = max_delay;    % CPP_len >= l_max-1
N_data = N-CPP_len;     % length of data symbols

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
gs = zeros(max_delay+1, N);
for i = 1:taps
    gs(delay_taps(i)+1, :) = gs(delay_taps(i)+1, :) + weighted(i, :);
end

% channel matrix form
H = Gen_channel_mtx(N, taps, chan_coef, delay_taps, Doppler_freq, c1);  % equation (24) in [R1]
% Observe the structure of H
% imagesc(abs(H))

% Loop Pre
gsConst = parallel.pool.Constant(gs);

ber_afdm = ones(size(SNR_dB));

%% loop begin
for iesn0 = 1:length(SNR_dB)
    sigma2 = sigma_2(iesn0);
    parErr = zeros(N_frame, 3);
    for iframe = 1:N_frame
        %% Random Data Generation
        info_bits = randi([0 1], N_data*log2(M_mod), 1);
        bits_reshape = reshape(info_bits, log2(M_mod), []).';
        symbol = bi2de(bits_reshape, 'left-msb');
        symbol_qam = qammod(symbol, M_mod, 'gray', 'UnitAveragePower', true);

        % Generate white noise
        w = sqrt(sigma2/2) * (randn(N, 1) + 1i*randn(N, 1));
        
        %% AFDM
        % Transmit
        s_afdm = AFDM_mod(symbol_qam, c1, c2);
        cpp_afdm = s_afdm(N_data-CPP_len:N_data-1).*exp(-1i*2*pi*c1*(N^2+2*N*(-CPP_len:-1).'));
        s_cpp_afdm = [cpp_afdm; s_afdm];
        % Through the channel
        gs_local = gsConst.Value;
        r_afdm = zeros(N,1);
        for l = (L_set+1)
            r_afdm(l:N) = r_afdm(l:N) + gs_local(l, l:N).' .* s_cpp_afdm(1:N-l+1);
        end
        r_afdm = r_afdm + w;
        % MMSE channel estimation
        x_est_afdm = H'/(H*H' + sigma2*eye(N)) * r_afdm;
        x_est_no_cpp_afdm = x_est_afdm(CPP_len+1:end);

        y_afdm = AFDM_demod(x_est_no_cpp_afdm, c1, c2);
        sym_det_afdm = qamdemod(y_afdm, M_mod, 'gray', 'UnitAveragePower', true);
        
        bits_det_afdm = de2bi(sym_det_afdm, log2(M_mod), 'left-msb').';
        bits_det_afdm = bits_det_afdm(:);
        err_afdm = sum(bits_det_afdm ~= info_bits);
        % fprintf('AFDM错误比特数 = %d\n', err_afdm);
        

        parErr(iframe, :) = [err_afdm, err_afdm, err_afdm];
    end  % parfor

    err_sum_AFDM = sum(parErr(:,1));

    ber_afdm(iesn0) = err_sum_AFDM / (log2(M_mod) * N_data * N_frame);

    fprintf('SNR=%2d dB done: afdm=%.3e\n', SNR_dB(iesn0), ber_afdm(iesn0));
end

% 绘图
figure;
semilogy(SNR_dB, ber_afdm, '-o', 'LineWidth', 1.2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
title(sprintf('N=%d taps=%d', N, taps));

toc
