% =========================================================================
% 复现 AFDM 论文 Fig. 9 的仿真结果
% Paper: Affine Frequency Division Multiplexing for Next Generation Wireless Communications
% 对应章节: VII. SIMULATION RESULTS
% =========================================================================

clc; clear; close all;

% --- 全局仿真参数 ---
N = 16;                 % 子载波数量 (Fig. 9)
mod_order = 2;          % BPSK
alpha_max = 1;          % 最大归一化多普勒频移 (Section VII)
SNR_dB_range = 0:2:20;  % SNR 范围

% 为了演示速度，减少了蒙特卡洛次数。论文中使用 10^6。
% 建议设置为至少 1000-5000 以观察大致趋势，设置为 10000+ 以获得平滑曲线。
num_frames = 1e2;

fprintf('开始仿真 Fig. 9a: 不同路径数下 AFDM 的性能...\n');
% === Fig. 9a 配置: AFDM, 不同路径数 P ===
% 2-path (l_max=1), 3-path (l_max=2), 4-path (l_max=3)
path_configs = [2, 1; 3, 2; 4, 3]; % [P, l_max]
ber_afdm_paths = zeros(size(path_configs, 1), length(SNR_dB_range));

for p_idx = 1:size(path_configs, 1)
    P = path_configs(p_idx, 1);
    l_max = path_configs(p_idx, 2);
    fprintf('正在仿真 P=%d, l_max=%d ...\n', P, l_max);

    ber_afdm_paths(p_idx, :) = run_simulation(N, P, l_max, alpha_max, ...
        SNR_dB_range, num_frames, "AFDM");
end

fprintf('开始仿真 Fig. 9b: 不同波形的性能对比...\n');
% === Fig. 9b 配置: 3-path channel, 对比不同波形 ===
P_comp = 3;
l_max_comp = 2;
waveforms = ["OFDM", "OCDM", "OTFS", "AFDM"];
ber_waveforms = zeros(length(waveforms), length(SNR_dB_range));

for w_idx = 1:length(waveforms)
    wf = waveforms(w_idx);
    fprintf('正在仿真波形: %s ...\n', wf);
    ber_waveforms(w_idx, :) = run_simulation(N, P_comp, l_max_comp, alpha_max, ...
        SNR_dB_range, num_frames, wf);
end

% === 绘图 ===
plot_results(SNR_dB_range, ber_afdm_paths, ber_waveforms, path_configs, waveforms);


function ber_curve = run_simulation(N, P, l_max, alpha_max, SNR_vec, num_frames, waveform)
% 主仿真循环
ber_curve = zeros(1, length(SNR_vec));

% 预计算发射符号的全排列 (用于 N=16 的 ML 检测)
% 注意: N=16, BPSK => 2^16 = 65536 种可能性。
% 预先生成所有可能的向量 x_all 用于暴力搜索 (为了代码简洁性)
% 在实际大 N 系统中应使用球形译码或近似算法。
num_candidates = 2^N;
% 使用十进制转二进制生成所有可能的 BPSK 向量 (-1, 1)
% 为了速度，这里仅在每次检测时针对 Heff 做处理，
% 或者我们可以使用简化的 Sphere Decoder，这里为了保证完全复现 ML 性能，
% 我们在 N=16 时使用蛮力搜索的优化版本。

% 生成 BPSK 候选集矩阵 (N x 2^N)
% x_candidates = 1 - 2*de2bi(0:num_candidates-1, N).';
% 上述生成可能内存过大或慢，这里我们采用分批或简单的逐帧检测。
% 鉴于 Matlab 矩阵运算优势，直接生成矩阵。
if N <= 16
    x_candidates = 1 - 2*dec2bin(0:num_candidates-1, N).' + 2*'0';
else
    error('N too large for brute force ML reproduction');
end

for snr_idx = 1:length(SNR_vec)
    snr = SNR_vec(snr_idx);
    N0 = 10^(-snr/10); % 假设信号能量归一化为 1

    total_errors = 0;
    total_bits = 0;

    parfor f = 1:num_frames % 使用并行计算加速
        % 1. 生成随机比特和符号
        bits = randi([0, 1], N, 1);
        x = 1 - 2*bits; % BPSK mapping: 0->1, 1->-1 (或者反过来，只要一致即可)

        % 2. 生成时域信道矩阵 H_time
        % 路径增益: Complex Gaussian, 0 mean, 1/P variance
        h = (randn(P, 1) + 1j*randn(P, 1)) / sqrt(2 * P);

        % 路径时延: 整数时延，分布在 [0, l_max]
        % 论文 implies distinct delays usually 0, 1, ... P-1
        delays = 0:(P-1);

        % 路径多普勒: Jakes Spectrum
        % alpha_i = alpha_max * cos(theta), theta ~ U[-pi, pi]
        thetas = (rand(P, 1) * 2 - 1) * pi;
        dopplers = alpha_max * cos(thetas);

        H_time = get_channel_time(N, h, delays, dopplers);

        % 3. 获取有效信道矩阵 H_eff 和 接收信号
        % y = H_eff * x + w

        U = eye(N); % 变换矩阵
        switch waveform
            case "AFDM"
                % Eq (47): c1 = (2*alpha_max + 1) / (2*N)
                c1 = (2*alpha_max + 1) / (2*N);
                % c2: 任意无理数或足够小。取黄金分割数相关小值。
                c2 = (sqrt(5)-1)/2 * (1/(20*N));
                A = get_daft_matrix(N, c1, c2);
                % AFDM: s = A' * x (Modulation), y_daft = A * r (Demodulation)
                % H_eff = A * H_time * A'
                U = A;

            case "OFDM"
                F = fft(eye(N)) / sqrt(N);
                % OFDM: s = F' * x, y_freq = F * r
                % H_eff = F * H_time * F'
                U = F;

            case "OCDM"
                c1 = 1/(2*N);
                c2 = 1/(2*N);
                A = get_daft_matrix(N, c1, c2);
                U = A;

            case "OTFS"
                % OTFS (SFFT based):
                % N = M * K. Let M=K=4 for N=16.
                M = 4; K = 4;
                % Effective channel in Delay-Doppler domain
                U = get_otfs_transform(M, K);
        end

        H_eff = U * H_time * U';

        % 4. 添加噪声
        % 信号经过幺正变换后能量不变
        w = (randn(N, 1) + 1j*randn(N, 1)) * sqrt(N0/2);
        y = H_eff * x + w;

        % 5. ML 检测
        % min || y - H_eff * x_cand ||^2
        % 展开: ||y||^2 + ||H x||^2 - 2 Re{y' H x}
        % 对于常模星座(BPSK/QPSK) 和 幺正变换近似，||H x||^2 项可能近似常数
        % 但在 LTV 信道下 H_eff 不是对角阵，我们做标准计算。

        % 矩阵化计算所有候选者的无噪声接收信号
        y_candidates = H_eff * x_candidates;

        % 计算欧氏距离
        dists = sum(abs(y - y_candidates).^2, 1);
        [~, min_idx] = min(dists);
        x_est = x_candidates(:, min_idx);

        % 统计误码
        % x_est 是 +1/-1, 转回 bits
        bits_est = (1 - x_est) / 2;
        num_err = sum(bits ~= bits_est);

        total_errors = total_errors + num_err;
        total_bits = total_bits + N;
    end

    ber_curve(snr_idx) = total_errors / total_bits;
    fprintf('  SNR = %d dB, BER = %.2e\n', snr, ber_curve(snr_idx));
end
end

% --- 辅助函数 ---

function H_time = get_channel_time(N, h, delays, dopplers)
% 根据公式 (24) 构建时域信道矩阵
% H = sum( h_i * Delta_fi * Pi^li )
% 忽略 CPP 矩阵，假设循环前缀/后缀完美处理，使得信道呈现循环卷积特性

H_time = zeros(N, N);

% 循环移位矩阵 Pi (Forward Cyclic Shift)
% Pi * x 将 x 向下移位 1 (x[n] -> x[n-1])
Pi = zeros(N, N);
Pi(1, N) = 1;
for k = 2:N
    Pi(k, k-1) = 1;
end

for i = 1:length(h)
    % Doppler Matrix Delta_fi = diag(exp(j*2pi*fi*n))
    % normalized doppler alpha = N * fi
    % term: exp(-j * 2pi * fi * n) in Eq (23) implies frequency shift
    % 注意：公式 (23) 是 exp(-j 2pi fi n)，公式 (24) Delta 定义是 diag(...)
    % 仔细对照 Eq (24) 和 (25):
    % Delta_fi = diag(exp(-j * 2pi * (alpha_i/N) * (0:N-1)))

    phase_shift = exp(-1j * 2 * pi * (dopplers(i)/N) * (0:N-1).');
    Delta_f = diag(phase_shift);

    % Matrix power for delay
    Pi_l = Pi^delays(i); % 实际上可以通过索引移位快速实现

    % path contribution
    H_time = H_time + h(i) * Delta_f * Pi_l;
end
end

function A = get_daft_matrix(N, c1, c2)
% 根据公式 (12) 生成 DAFT 变换矩阵 (Receiver side transform)
% Sm = (1/sqrt(N)) * sum_n exp(-j2pi(c2*m^2 + 1/N*m*n + c1*n^2)) * sn
% 矩阵元素 A(m,n) 对应 m (0..N-1) 行, n (0..N-1) 列

[n_grid, m_grid] = meshgrid(0:N-1, 0:N-1);

% Eq (12) 指数项
phase = 2 * pi * (c2 * m_grid.^2 + (1/N) * m_grid .* n_grid + c1 * n_grid.^2);
A = (1/sqrt(N)) * exp(-1j * phase);
end

function U_otfs = get_otfs_transform(M, K)
% 构建 OTFS 的等效变换矩阵 U (Time -> Delay-Doppler)
% 使得 y_dd = U * y_time
% 实际上，OTFS 检测通常在 DD 域。
% 关系: x_time = Heisenberg(x_dd)
% y_time = H_time * x_time
% y_dd = Wigner(y_time)
% 所以 H_eff_dd = Wigner * H_time * Heisenberg
% 此处 U = Wigner 变换矩阵。

% ISFFT (DD to TF): F_M * X_dd * F_K' (通常定义)
% Heisenberg (TF to Time): 相当于对 TF 数据做 IFFT (沿频率轴) ?
% 为了简单且标准化，使用 Kronecker 积定义 OTFS 调制矩阵 (DD向量 -> 时域向量)
% 引用常用模型: x_time = (kron(F_M', I_K)) * x_dd (假设 x_dd 按列堆叠)
% 或者依据 SFFT 定义。

% 这里我们使用幺正变换定义：
% U_otfs 将 时域 映射到 DD 域。
% Isff_matrix: DD -> Time
% Isff = kron(F_M', eye(K)); % 这是一个简化的变换，实际上 OTFS 还要更复杂

% 按照论文描述：SFFT 实现。
% 发送端 (DD -> Time): x_dd -> F_M(row) -> F_K'(col) -> Heisenberg
% 简易实现：
% 1. ISFFT: DD -> TF. X_tf(n,m) = sum X_dd ...
% 2. Heisenberg: TF -> Time.

% 构建显式矩阵
N = M * K;
I_mat = eye(N);

% 构造解调矩阵 U (Time -> DD)
% 逐列处理单位阵，看它们变换后的结果，组合成矩阵
U_otfs = zeros(N, N);

F_M = fft(eye(M)) / sqrt(M);
F_K = fft(eye(K)) / sqrt(K);

for col = 1:N
    y_time_vec = I_mat(:, col);
    y_time_mat = reshape(y_time_vec, K, M); % Assuming column-major time framing

    % Wigner Transform (Time -> TF)
    % Standard OFDM demod (Time -> Freq) per symbol
    Y_tf = fft(y_time_mat, K, 1) / sqrt(K); % K-point FFT on columns

    % SFFT (TF -> DD)
    % DD = SFFT(TF) = F_K * TF * F_M' ?
    % 论文引用 [28] SFFT: F_M^H \otimes F_N (注意符号)
    % 通常 SFFT 是：沿多普勒轴(M)做 FFT，沿时延轴(K)做 IFFT。
    % 但 OTFS 中，接收端是从 TF 到 DD。
    % TF(k, m) -> DD(l, k_doppler)
    % 这是一个辛傅里叶变换。
    % 简单起见，既然 AFDM 也是幺正变换，OTFS 也是。
    % 使用标准 OTFS 矩阵定义: U = (F_M x I_K) * P ?
    % 我们采用 Zak 变换形式的数值计算:

    % Wigner / Demod to DD:
    % 1. Reshape to K x M
    % 2. FFT along K (Delay axis processing from Time) -> Y_tf
    % 3. IFFT along M (Doppler axis recovery)
    Y_dd = ifft(Y_tf, M, 2) * sqrt(M);

    U_otfs(:, col) = Y_dd(:);
end
end

function plot_results(snr, ber_afdm, ber_comp, path_configs, waveforms)
% 绘制 Fig 9a
figure;
semilogy(snr, ber_afdm(1, :), 'r-o', 'LineWidth', 1.5, 'DisplayName', 'AFDM (2 paths)');
hold on;
semilogy(snr, ber_afdm(2, :), 'b-s', 'LineWidth', 1.5, 'DisplayName', 'AFDM (3 paths)');
semilogy(snr, ber_afdm(3, :), 'g-^', 'LineWidth', 1.5, 'DisplayName', 'AFDM (4 paths)');

% 添加参考斜率 (如 s ~ SNR^-P)
% 仅示意
ref_snr = snr(end-2:end);
ref_val = ber_afdm(2, end-2:end);
semilogy(ref_snr, ref_val(1) * (10.^(-(0:2)*3/10)), 'k--', 'DisplayName', 'Slope Order 3');

grid on;
xlabel('SNR (dB)'); ylabel('BER');
title('Fig. 9a: AFDM BER with different number of paths');
legend('Location', 'southwest');

% 绘制 Fig 9b
figure;
markers = {'d--', 'v--', 'o-', 's-'};
colors = {'r', 'm', 'b', 'g'};
for i = 1:length(waveforms)
    semilogy(snr, ber_comp(i, :), [colors{i} markers{i}], 'LineWidth', 1.5, ...
        'DisplayName', waveforms(i));
    hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('BER');
title('Fig. 9b: Comparison of waveforms (3-path LTV channel)');
legend('Location', 'southwest');
end