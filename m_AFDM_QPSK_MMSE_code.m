%% FOR Fig. 10. BER and spectral efficiency performance of OFDM, OCDM, OTFS and AFDM using MMSE detection.

clear;
clc;

% rng(7)
tic
% System parameters
M_mod = 4;      % size of QAM constellation
N = 256;        % number of symbols(subcarriers)
B = 10e6;

N_frame = 10;    % number of simulation frames

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


% k = log2(M_mod);
% Rc = 1/2;                % 卷积码 rate
% EbN0_dB = SNR_dB;       % 把原来的 SNR_dB 解释为 Eb/N0（推荐这样比较）
% EbN0 = 10.^(EbN0_dB/10);
% 
% EsN0 = EbN0 .* (k * Rc);   % Es/N0 对应值
% sigma_2 = 1 ./ EsN0;        % Es=1 时，sigma2 = N0, noise power per complex sample = N0





% Generate synthetic delay-Doppler channel %% 生成合成延迟-多普勒信道
k_max = 2;    %  maximum Doppler shift is αmax = 2
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
fprintf("Doppler_freq : %.2fkHz.\n", Doppler_freq);

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

trellis = poly2trellis(7, [171 133]);

ber_uncoded = zeros(size(SNR_dB));
ber_coded   = zeros(size(SNR_dB));


%% loop begin
for iesn0 = 1:length(SNR_dB)
    sigma2 = sigma_2(iesn0);
    sigma2c = sigma_2_code(iesn0);
    parErr = zeros(N_frame, 2);
    for iframe = 1:N_frame
        %% Random Data Generation
        info_bits = randi([0 1], N_data*log2(M_mod), 1);
        % Generate white noise
        w = sqrt(sigma2/2) * (randn(N, 1) + 1i*randn(N, 1));
        w_c = sqrt(sigma2c/2) * (randn(N, 1) + 1i*randn(N, 1));
        

        %% AFDM uncode
        % uncoded: data Modulation
        bits_reshape = reshape(info_bits, log2(M_mod), []).';
        sym_uncoded = bi2de(bits_reshape, 'left-msb');
        x_qam_uncoded = qammod(sym_uncoded, M_mod, 'gray', 'UnitAveragePower', true);
        % x_qam_uncoded = qammod(sym_uncoded, M_mod, 'gray');


        % uncoded: data Transmit by AFDM
        s_afdm_unc = AFDM_mod(x_qam_uncoded, c1, c2);
        cpp_afdm_unc = s_afdm_unc(N_data-CPP_len:N_data-1).*exp(-1i*2*pi*c1*(N^2+2*N*(-CPP_len:-1).'));
        s_cpp_afdm_unc = [cpp_afdm_unc; s_afdm_unc];

        gs_local = gsConst.Value;
        r_afdm_unc = zeros(N,1);
        for l = (L_set+1)
            r_afdm_unc(l:N) = r_afdm_unc(l:N) + gs_local(l, l:N).' .* s_cpp_afdm_unc(1:N-l+1);
        end
        r_afdm_unc = r_afdm_unc + w;

        % MMSE
        x_est_afdm_unc = H'/(H*H' + sigma2*eye(N)) * r_afdm_unc;
        x_est_no_cpp_afdm_unc = x_est_afdm_unc(CPP_len+1:end);

        y_afdm_unc = AFDM_demod(x_est_no_cpp_afdm_unc, c1, c2);
        sym_det_unc = qamdemod(y_afdm_unc, M_mod, 'gray', 'UnitAveragePower', true);
        % sym_det_unc = qamdemod(y_afdm_unc, M_mod, 'gray');

        bits_det_unc = de2bi(sym_det_unc, log2(M_mod), 'left-msb').';
        bits_det_unc = bits_det_unc(:);
        err_unc = sum(bits_det_unc ~= info_bits);
        % fprintf('未编码的错误比特数 = %d\n', err_unc);

        %% AFDM code
        info_bits = randi([0 1], N_data*log2(M_mod)/2, 1);
        % 卷积编码
        coded_bits = convenc(info_bits, trellis); 
        % 交织
        matrix = reshape(coded_bits, log2(M_mod), []);
        intlvddata = matintrlv(matrix, 2, log2(M_mod) / 2);
        % QAM 调制
        sym_coded =  bi2de(intlvddata','left-msb');
        % x_qam_coded = qammod(sym_coded, M_mod, 'gray');
        x_qam_coded = qammod(sym_coded, M_mod, 'gray', 'UnitAveragePower', true);

% % 在发送端直接按 k 分组并 bi2de (不交织)
% sym_coded = bi2de(reshape(coded_bits, log2(M_mod), []).','left-msb');
% x_qam_coded = qammod(sym_coded, M_mod, 'gray', 'UnitAveragePower', true);

        % AFDM 调制
        s_afdm_coded = AFDM_mod(x_qam_coded, c1, c2);
        cpp_afdm_coded = s_afdm_coded(N_data-CPP_len:N_data-1).*exp(-1i*2*pi*c1*(N^2+2*N*(-CPP_len:-1).'));
        s_cpp_afdm_coded = [cpp_afdm_coded; s_afdm_coded];
        % 信道衰落
        r_afdm_coded = zeros(N,1);
        for l = (L_set+1)
            r_afdm_coded(l:N) = r_afdm_coded(l:N) + gs_local(l, l:N).' .* s_cpp_afdm_coded(1:N-l+1);
        end
        % 信道噪声
        r_afdm_coded = r_afdm_coded + w_c;
        % MMSE 信道检测
        x_est_afdm_coded = H'/(H*H' + sigma2c*eye(N)) * r_afdm_coded;   % MMSE
        x_est_no_cpp_afdm_coded = x_est_afdm_coded(CPP_len+1:end);
        % AFDM 解调
        y_afdm_coded = AFDM_demod(x_est_no_cpp_afdm_coded, c1, c2);
        % QAM 解调
        % sym_det_coded = qamdemod(y_afdm_coded, M_mod, 'gray');
        sym_det_coded = qamdemod(y_afdm_coded, M_mod, 'gray', 'UnitAveragePower', true);
        bits_det_coded = de2bi(sym_det_coded, log2(M_mod), 'left-msb').';
         % 反交织
        deinterleaved_bits = matdeintrlv(bits_det_coded, 2, log2(M_mod) / 2);
        deinterleaved_bits = deinterleaved_bits(:);
        % 卷积译码
        % decoded_bits = vitdec(deinterleaved_bits', trellis, 25, 'trunc', 'unquant')';
        decoded_bits = vitdec(deinterleaved_bits', trellis, 25, 'trunc', 'hard')';

% % 在接收端，直接 de2bi 后得到 bits_rx 列向量，将其传给 vitdec
% bits_det_coded = de2bi(sym_det_coded, log2(M_mod), 'left-msb').'; bits_det_coded = bits_det_coded(:);
% decoded_bits = vitdec(bits_det_coded', trellis, 35, 'trunc', 'hard')';   % 注意转置

        err_coded = sum(decoded_bits ~= info_bits);

        err_pair = [err_unc, err_coded];
        parErr(iframe, :) = err_pair;
    end  % parfor

    err_sum_uncoded = sum(parErr(:,1));
    err_sum_coded = sum(parErr(:,2));

    ber_uncoded(iesn0) = err_sum_uncoded / (log2(M_mod) * N_data * N_frame);
    ber_coded(iesn0)   = err_sum_coded   / (log2(M_mod) * N_data / 2 * N_frame); 

    fprintf('SNR=%2d dB done: uncoded=%.3e, coded=%.3e\n', SNR_dB(iesn0), ber_uncoded(iesn0), ber_coded(iesn0));
end

% ---- 绘图比较（按 information-bit 的 BER 比较） ----
figure;
semilogy(SNR_dB, ber_uncoded, '-o', 'LineWidth', 1.2); hold on;
semilogy(SNR_dB, ber_coded, '-s', 'LineWidth', 1.2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
legend('Uncoded','Coded (conv + interleaver)');
title(sprintf('N=%d taps=%d', N, taps));

toc
