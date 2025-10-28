%% FOR Fig. 10. BER and spectral efficiency performance of OFDM, OCDM, OTFS and AFDM using MMSE detection.

clear; clc;
% close all;
rng(7)
tic
%% System parameters %%
M_mod = 4;      % size of QAM constellation
N = 256;        % number of symbols(subcarriers)
B = 10e6;

car_fre = 4e9;  % carrier frequency
delta_f = 15e3;  % symbol spacing    符号间距
% delta_f = B/N;  % symbol spacing    符号间距
T = 1/delta_f;  % symbol duration   符号持续时间

eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));   % average power per symbol
SNR_dB = 0:2:20;    % set SNR here
SNR = 10.^(SNR_dB/10);
% sigma_2 = 1 ./ SNR;   % noise power
sigma_2 = (abs(eng_sqrt)^2)./SNR;   % noise power

N_frame = 1000;    % number of simulation frames

fprintf("Number of subcarriers : %d.\n", N);
fprintf("Total bandwidth       : %.2fMHz.\n", B/1e6);
fprintf("Symbol spacing        : %.2fkHz.\n", delta_f/1E3);
fprintf("Symbol duration       : %.2fus.\n", T*1e6);

%% Generate synthetic delay-Doppler channel %% 生成合成延迟-多普勒信道

%%% The maximum Doppler shift is αmax = 2,
%%% which corresponds to a speed of 540 km/h,
%%% and the Doppler shift of each path
%%% is generated using Jakes Doppler spectrum
k_max = 0;                  % maximum normalized Doppler index

%%% We consider a 3-path channel.
taps  = 1;        % number of paths

%%% The maximum delay spread is set to be lmax = 2.
l_max = 1;    % maximum normalized delay index

chan_coef = 1/sqrt(2).*(randn(1,taps)+1i.*randn(1,taps));   % follows Rayleigh distribution
fprintf("P = %d. Channel Power = %.2f.\n", taps, sum(abs(chan_coef).^2))

%%% integer delay shifts: random delays in range [0,l_max-1]
delay_taps = randi(l_max, [1,taps]) - 1;
% delay_taps = sort(delay_taps-min(delay_taps));
fprintf("delay_taps:"); disp(delay_taps);

%%% fractional Doppler shifts: uniformly distributed Doppler shifts in range [-k_max,k_max]
% Doppler_taps = k_max*(2*rand(1,taps)-1);
% % Doppler_taps = round(Doppler_taps);     % cast to integer Doppler shifts
% Doppler_freq = Doppler_taps/(N*T);      % f=k/(NT),f:Doppler shifts(Hz),k:normalized Doppler shifts
%%% the Doppler shift of each path is generated using Jakes Doppler spectrum
fD_max = k_max / (N*T);  % 最大物理多普勒频移
u = rand(1, taps);
Doppler_freq = fD_max * sin(pi * (u - 0.5));   % 服从近似Jakes分布
Doppler_taps = Doppler_freq * N*T;
fprintf("Doppler_freq : %.2fkHz.\n", Doppler_freq); 


%% AFDM parameters %%
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

%% Generate channel matrix %%
% discrete-time channel 离散时间信道
L_set = unique(delay_taps);
gs=zeros(max_delay+1,N);
for q=0:N-1
    for i=1:taps
        h_i=chan_coef(i);   % the complex gain
        l_i=delay_taps(i);  % the integer delay associated with the i-th path,
        f_i=Doppler_freq(i);% Doppler shift (in digital frequencies)
        % Dirac delta function 在零点以外的所有位置值为零，而在整个定义域上的积分值为1
        gs(l_i+1,q+1)=gs(l_i+1,q+1)+h_i*exp(-1i*2*pi*f_i*q);  % equation (23) in [R1]
    end
end


% channel matrix form
H = Gen_channel_mtx(N, taps, chan_coef, delay_taps, Doppler_freq, c1);  % equation (24) in [R1]
% Observe the structure of H
% imagesc(abs(H))

H_FIR = zeros(N,1);
for idx = 1:taps
    d = delay_taps(idx);   % integer delay
    H_FIR(d+1) = H_FIR(d+1) + chan_coef(idx);  % 注意下标 (1-based)
end

%% Start Loop
ber_AFDM  = zeros(size(SNR_dB));
ber_OFDM  = zeros(size(SNR_dB));
ber_OCDM  = zeros(size(SNR_dB));

for iesn0 = 1:length(SNR_dB)

    sigma2 = sigma_2(iesn0);
    err_count_AFDM = zeros(N_frame,1);
    err_count_OFDM = zeros(N_frame,1);
    err_count_OCDM = zeros(N_frame,1);

    parfor iframe = 1:N_frame
        %% Tx data generation %%
        x = randi([0, M_mod-1], N_data, 1);     % generate random bits
        x_qam = qammod(x, M_mod, 'gray', 'UnitAveragePower', true);
        w = sqrt(sigma2/2) * (randn(N, 1) + 1i*randn(N, 1));

        %% ========== AFDM chain ==========
        %%% AFDM modulation
        s_afdm = AFDM_mod(x_qam, c1, c2);
        %%% generate CPP
        cpp_afdm = s_afdm(N_data-CPP_len:N_data-1).*exp(-1i*2*pi*c1*(N^2+2*N*(-CPP_len:-1).'));
        %%% Insert CPP
        s_cpp_afdm = [cpp_afdm; s_afdm];
        %%% Through dela y-Doppler channel
        r_afdm = zeros(N,1);
        for l = (L_set+1)
            r_afdm(l:N) = r_afdm(l:N) + gs(l, l:N).' .* s_cpp_afdm(1:N-l+1);
        end
        %%% Add AWGN
        r_afdm = r_afdm + w;
        %%% MMSE equalization, ideal channel estimation
        x_est_afdm = H'/(H*H'+sigma2*eye(N))*r_afdm;
        %%% discard CPP
        x_est_no_cpp_afdm = x_est_afdm(CPP_len+1:end);
        %%% AFDM demodulation
        y_afdm = AFDM_demod(x_est_no_cpp_afdm, c1, c2);
        x_est_bit_afdm = qamdemod(y_afdm, M_mod, 'gray');
        %%% Error count %%
        err_count_AFDM(iframe) = sum(x_est_bit_afdm ~= x);

        %% ========== OFDM chain ==========
        S_ofdm = ifft(x_qam, N_data_ofdm) * sqrt(N_data_ofdm/N);
        cpp_ofdm = S_ofdm(end-CP_len+1:end);
        scpp_ofdm = [cpp_ofdm; S_ofdm];

        r_ofdm = zeros(N,1);
        for n = 1:N
            for i = 1:taps
                delay = delay_taps(i);
                doppler = Doppler_freq(i);
                if n > delay
                    r_ofdm(n) = r_ofdm(n) + chan_coef(i)*exp(1j*2*pi*doppler*n*T)*scpp_ofdm(n-delay);
                end
            end
        end
        r_ofdm = r_ofdm + w;

        % % OFDM modulation: IFFT of data (length N_data)
        % S_ofdm = ifft(x_qam, N_data);
        % % OFDM cyclic prefix: last CPP_len samples copied
        % cpp_ofdm = S_ofdm(end-CPP_len+1:end);
        % scpp_ofdm = [cpp_ofdm; S_ofdm];
        % r_ofdm = zeros(N,1);
        % for l=(L_set+1) 
        %     r_ofdm(l:N) = r_ofdm(l:N) + gs(l, l:N).' .* scpp_ofdm(1:N-l+1);
        % end
        % r_ofdm = r_ofdm + w;
        r_ofdm_no_cp = r_ofdm(CPP_len+1 : CPP_len+N_data);
        % Frequency-domain channel on N_data points
        H_f = fft(H_FIR, N_data);    % length N_data
        % FFT the received OFDM time-block
        Y_k = fft(r_ofdm_no_cp, N_data);
        % MMSE per-subcarrier equalizer (scalar)
        Xhat_k = (conj(H_f) ./ (abs(H_f).^2 + sigma2)) .* Y_k;
        y_ofdm = Xhat_k;   
        x_est_bit_ofdm = qamdemod(y_ofdm, M_mod, 'gray');
        err_count_OFDM(iframe) = sum(x_est_bit_ofdm ~= x);

        % % % % % OFDM modulation: IFFT of data (length N_data)
        % % % % S_ofdm = ifft(x_qam, N_data);
        % % % % % OFDM cyclic prefix: last CPP_len samples copied
        % % % % cpp_ofdm = S_ofdm(end-CPP_len+1:end);
        % % % % scpp_ofdm = [cpp_ofdm; S_ofdm];
        % % % % r_ofdm = zeros(N,1);
        % % % % for l=(L_set+1)
        % % % %     r_ofdm(l:N) = r_ofdm(l:N) + gs(l, l:N).' .* scpp_ofdm(1:N-l+1);
        % % % % end
        % % % % r_ofdm = r_ofdm + w;
        % % % % 
        % % % % r_ofdm_no_cp = r_ofdm(CPP_len+1 : CPP_len+N_data);
        % % % % 
        % % % % % Frequency-domain channel on N_data points
        % % % % H_f = fft(H_FIR, N_data);    % length N_data
        % % % % 
        % % % % % FFT the received OFDM time-block
        % % % % Y_k = fft(r_ofdm_no_cp, N_data);
        % % % % 
        % % % % % MMSE per-subcarrier equalizer (scalar)
        % % % % % 注意 sigma2 已在外层定义为噪声方差
        % % % % % MMSE: Xhat_k = conj(H_k) ./ (|H_k|^2 + sigma2) .* Y_k
        % % % % Xhat_k = (conj(H_f) ./ (abs(H_f).^2 + sigma2)) .* Y_k;
        % % % % 
        % % % % % recover QAM symbols and demodulate
        % % % % y_ofdm = Xhat_k;   % frequency-domain estimated symbols
        % % % % x_est_bit_ofdm = qamdemod(y_ofdm, M_mod, 'gray');
        % % % % % x_est_ofdm = H'/(H*H'+sigma2*eye(N))*r_ofdm;
        % % % % % x_est_no_cpp_ofdm = x_est_ofdm(CPP_len+1:end);
        % % % % % % OFDM demod: FFT
        % % % % % y_ofdm = fft(x_est_no_cpp_ofdm, N_data);
        % % % % % x_est_bit_ofdm = qamdemod(y_ofdm, M_mod, 'gray');
        % % % % err_count_OFDM(iframe) = sum(x_est_bit_ofdm ~= x);

        %% ========== OCDM chain (DFnT-based) ==========
        s_ocdm = OCDM_mod(x_qam);    % returns length N_data
        % For OCDM, we'll use conventional CP (copy last CPP_len)
        cpp_ocdm = s_ocdm(end-CPP_len+1:end);
        scpp_ocdm = [cpp_ocdm; s_ocdm];
        r_ocdm = zeros(N,1);
        for l=(L_set+1)
            r_ocdm(l:N) = r_ocdm(l:N) + gs(l, l:N).' .* scpp_ocdm(1:N-l+1);
        end
        r_ocdm = r_ocdm + w;
        x_est_ocdm = H'/(H*H'+sigma2*eye(N))*r_ocdm;
        x_est_no_cpp_ocdm = x_est_ocdm(CPP_len+1:end);
        % OCDM demod: inverse DFnT
        y_ocdm = OCDM_demod(x_est_no_cpp_ocdm);
        x_est_bit_ocdm = qamdemod(y_ocdm, M_mod, 'gray');
        err_count_OCDM(iframe) = sum(x_est_bit_ocdm ~= x);
        
    end
    ber_AFDM(iesn0) = sum(err_count_AFDM)/N_data/N_frame;
    ber_OFDM(iesn0) = sum(err_count_OFDM)/N_data/N_frame;
    ber_OCDM(iesn0) = sum(err_count_OCDM)/N_data/N_frame;

    fprintf('SNR=%2d dB done: AFDM=%.3e OFDM=%.3e OCDM=%.3e\n', SNR_dB(iesn0), ber_AFDM(iesn0), ber_OFDM(iesn0), ber_OCDM(iesn0));

end

%% Plot bit error rate %%
figure;
semilogy(SNR_dB, ber_AFDM, '-o', 'LineWidth', 1.2); hold on;
semilogy(SNR_dB, ber_OFDM, '-s', 'LineWidth', 1.2);
semilogy(SNR_dB, ber_OCDM, '-^', 'LineWidth', 1.2);
grid on;
xlabel('SNR (dB)');
ylabel('Symbol error rate (per symbol)');
title(sprintf("N=%d P=%d K=%d L=%d", N, taps, k_max, l_max));
legend('AFDM','OFDM','OCDM','Location','southwest');

toc

%% ========== Inline helper functions ==========

function X = OCDM_mod(x)
% OCDM modulation via a DFnT-like transform (length = length(x))
Nloc = length(x);
n = (0:Nloc-1).';
chirp = exp(1i*pi*(n.^2)/Nloc);
X = chirp .* fft(chirp .* x);
% normalize to unit average power similar to IFFT scaling (optional)
X = X / sqrt(mean(abs(X).^2)) * sqrt(mean(abs(x).^2));
end

function x = OCDM_demod(X)
Nloc = length(X);
n = (0:Nloc-1).';
chirp = exp(1i*pi*(n.^2)/Nloc);
x = conj(chirp) .* ifft(conj(chirp) .* X);
% keep power consistent
x = x / sqrt(mean(abs(x).^2)) * sqrt(mean(abs(X).^2));
end