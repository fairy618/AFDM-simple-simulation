%% FOR Fig. 10. BER and spectral efficiency performance of OFDM, OCDM, OTFS and AFDM using MMSE detection.

clear; clc;
% close all;
rng(1)

%% System parameters %%
M_mod = 4;      % size of QAM constellation
N = 256;        % number of symbols(subcarriers)

car_fre = 4e9;  % carrier frequency
delta_f = 1e3;  % symbol spacing    符号间距
T = 1/delta_f;  % symbol duration   符号持续时间

eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));   % average power per symbol
SNR_dB = 0:2:20;    % set SNR here
SNR = 10.^(SNR_dB/10);
sigma_2 = 1 ./ SNR;   % noise power

N_frame = 10000;    % number of simulation frames

%% Generate synthetic delay-Doppler channel %% 生成合成延迟-多普勒信道

%%% The maximum Doppler shift is αmax = 2,
%%% which corresponds to a speed of 540 km/h,
%%% and the Doppler shift of each path
%%% is generated using Jakes Doppler spectrum
k_max = 2;                  % maximum normalized Doppler index

%%% We consider a 3-path channel.
taps  = 3;        % number of paths

%%% The maximum delay spread is set to be lmax = 2.
l_max = 2;    % maximum normalized delay index

% chan_coef = 1/sqrt(2).*(randn(1,taps)+1i.*randn(1,taps));   % follows Rayleigh distribution
chan_coef = (randn(1,taps)+1i*randn(1,taps));
chan_coef = chan_coef / norm(chan_coef);   % 归一化总功率为1

fprintf("P = %d. Channel Power = %.2f.\n", taps, sum(abs(chan_coef).^2))

%%% integer delay shifts: random delays in range [0,l_max-1]
delay_taps = randi(l_max, [1,taps]) - 1;
% delay_taps = sort(delay_taps-min(delay_taps));

%%% fractional Doppler shifts: uniformly distributed Doppler shifts in range [-k_max,k_max]
% Doppler_taps = k_max*(2*rand(1,taps)-1);
% % Doppler_taps = round(Doppler_taps);     % cast to integer Doppler shifts
% Doppler_freq = Doppler_taps/(N*T);      % f=k/(NT),f:Doppler shifts(Hz),k:normalized Doppler shifts
%%% the Doppler shift of each path is generated using Jakes Doppler spectrum
fD_max = k_max / (N*T);  % 最大物理多普勒频移
u = rand(1, taps);
Doppler_freq = fD_max * sin(pi * (u - 0.5));   % 服从近似Jakes分布
Doppler_taps = Doppler_freq * N*T;

% fprintf("P=%d, Channel Power=%.2f\n", taps, sum(abs(chan_coef).^2));

%% AFDM parameters %%
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


figure();

ber = zeros(size(SNR_dB));
for iesn0 = 1:length(SNR_dB)

    sigma2 = sigma_2(iesn0);

    err_count = zeros(size(N_frame));
    for iframe = 1:N_frame
        %% Tx data generation %%
        x = randi([0, M_mod-1], N_data, 1);     % generate random bits

        % x_qam = qammod(x, M_mod, 'gray');   % QAM modulation
        x_qam = qammod(x, M_mod, 'gray', 'UnitAveragePower', true);

        s = AFDM_mod(x_qam, c1, c2);    % AFDM modulation

        cpp = s(N_data-CPP_len:N_data-1).*exp(-1i*2*pi*c1*(N^2+2*N*(-CPP_len:-1).'));     % generate CPP
        s_cpp = [cpp; s];   % Insert CPP

        %% Through delay-Doppler channel %%
        r=zeros(N,1);
        for q=1:N
            for l=(L_set+1)
                if(q>=l)
                    r(q)=r(q)+gs(l,q)*s_cpp(q-l+1);  %equation (22) in [R1]
                end
            end
        end
        w = sqrt(sigma2/2) * (randn(size(s_cpp)) + 1i*randn(size(s_cpp)));    % add Gaussian noise
        r=r+w;
        % r=H*s_cpp+w;  % or simply do this

        %% Rx detection %%
        x_est = H'/(H*H'+sigma2*eye(N))*r;  % MMSE equalization, ideal channel estimation

        x_est_no_cpp = x_est(CPP_len+1:end);  % discard CPP

        y = AFDM_demod(x_est_no_cpp, c1, c2);  % AFDM demodulation

        x_est_bit = qamdemod(y, M_mod, 'gray');  % QAM demodulation

        %% Error count %%
        err_count(iframe) = sum(x_est_bit ~= x);    % calculate error bits
    end
    ber(iesn0) = sum(err_count)/length(x)/N_frame;  % calculate bit error rate
end

disp(ber)

%% Plot bit error rate %%
semilogy(SNR_dB, ber)

xlabel('SNR(dB)')
ylabel('BER')
title('BER of AFDM systems')
grid on


