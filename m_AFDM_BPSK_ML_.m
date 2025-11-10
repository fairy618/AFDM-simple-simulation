% Fig 9a

clear; clc;
% close all;
rng(1)

%% System parameters %%
% The carrier frequency is 4 GHz.
fc = 4e9;   % carrier frequency
c = 3e8;
lambda = c/fc;
% The duration between two successive delay taps is approximately 41.6 μs
T_tap = 41.6e-6;
fs = 1/T_tap;

% N = 16 and BPSK
M_mod = 2;
N = 16;

T_sym = N * T_tap;
delta_f = 1/T_sym;


eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));   % average power per symbol
SNR_dB = 0:2:20;    % set SNR here
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;   % noise power

N_frame = 1e5;    % number of simulation frames

%% Generate synthetic delay-Doppler channel %% 生成合成延迟-多普勒信道

figure();
for num_of_Path = 2:4

    k_max = 1;                  % maximum normalized Doppler index
    taps  = num_of_Path;        % number of paths
    l_max = num_of_Path - 1;    % maximum normalized delay index

    ber = zeros(size(SNR_dB));
    for iesn0 = 1:length(SNR_dB)

        sigma2 = sigma_2(iesn0);

        err_count = zeros(size(N_frame));
        parfor iframe = 1:N_frame

            %% Generate synthetic delay-Doppler channel %% 生成合成延迟-多普勒信道
            % follows Rayleigh distribution
            chan_coef = 1/sqrt(2*taps).*(randn(1,taps)+1i.*randn(1,taps));  % 总能量固定
            % chan_coef = 1/sqrt(2).*(randn(1,taps)+1i.*randn(1,taps));       % 每条路径独立随机

            delay_taps = randi(l_max, [1,taps]) - 1;
            
            fD_max = k_max / T_sym;  % 最大物理多普勒频移
            theta = (rand(1,taps)*2 - 1) * pi;   % uniform in [-pi, pi]
            Doppler_freq = fD_max * cos(theta);         % Hz

            Doppler_taps = Doppler_freq * T_sym;


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
            qq = 0:N-1; % 所有频率索引
            phase = exp(-1i*2*pi*(Doppler_freq(:) * qq));   % taps × N
            weighted = chan_coef(:) .* phase;  % taps × N
            gs=zeros(max_delay+1,N);
            for i = 1:taps
                gs(delay_taps(i)+1, :) = gs(delay_taps(i)+1, :) + weighted(i, :);
            end

            % for q=0:N-1
            %     for i=1:taps
            %         h_i=chan_coef(i);   % the complex gain
            %         l_i=delay_taps(i);  % the integer delay associated with the i-th path,
            %         f_i=Doppler_freq(i);% Doppler shift (in digital frequencies)
            %         % Dirac delta function 在零点以外的所有位置值为零，而在整个定义域上的积分值为1
            %         gs(l_i+1,q+1)=gs(l_i+1,q+1)+h_i*exp(-1i*2*pi*f_i*q);  % equation (23) in [R1]
            %     end
            % end

            % channel matrix form
            H = Gen_channel_mtx(N, taps, chan_coef, delay_taps, Doppler_freq, c1);  % equation (24) in [R1]
            % Observe the structure of H
            % imagesc(abs(H))


            %% Tx data generation %%
            x = randi([0, M_mod-1], N_data, 1);     % generate random bits

            x_qam = qammod(x, M_mod, 'gray');   % QAM modulation

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
        ber(iesn0) = sum(err_count)/N/N_frame;  % calculate bit error rate
    end

    % Plot bit error rate %%
    fprintf("P=%d, ", num_of_Path);
    disp(ber)
    semilogy(SNR_dB, ber)
    hold on

end
hold off

legend("P=2", "P=3", "P=4")

xlabel('SNR(dB)')
ylabel('BER')
title('BER of AFDM systems')
grid on


