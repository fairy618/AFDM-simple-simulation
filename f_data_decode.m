% 好的，根据你刚刚生成的f_data_code函数，还需要一个f_data_decode函数。
% 我认为f_data_decode函数需要的输入参数有这些：
% 1、rec_data： 接收到的调制信号（列向量）
% 2、M_mod是QAM调制阶数，例如4, 16, 64等（暂时不考虑非标准的，认为我输入的就是2、4、16、64等等）
% 3、trellis是我生成的卷积码的trellis结构。
% 我需要你做的就是，补充函数，包括：
% 1、将数据反交织
% 2、将反交织后的的数据重新整形为列向量
% 3、将整形后的数据根据trellis进行解码，恢复原始数据

% 我需要实现一个matlab的子函数：function y = f_data_code(data, M_mod, trellis)
% 其中data是需要编码的二进制数据（列向量）；
% M_mod是QAM调制阶数，例如4, 16, 64等（暂时不考虑非标准的，认为我输入的就是2、4、16、64等等）
% trellis是我生成的卷积码的trellis结构。
% 我需要你做的就是，补充函数，包括：
% 1、将数据卷积编码
% 2、将编码的数据重新整形
% 3、将整形后的数据进行交织

function y = f_data_decode(rec_data, M_mod, trellis)
% f_data_decode : 反交织 + 整形 + 卷积译码（硬判决）
% 输入：
%   rec_data : 接收到的比特列向量（0/1，已由 QAM 解调并量化）
%   M_mod    : QAM阶数（2,4,16,64,...）
%   trellis  : poly2trellis 返回的 trellis 结构
% 输出：
%   y : 解码恢复的原始比特列向量

rec_data = rec_data(:);   % 确保列向量
QAMbit = log2(M_mod);
if QAMbit ~= round(QAMbit)
    error('M_mod must be power of 2');
end

% 交织参数自动选择（与编码端 f_data_code 保持一致的逻辑）
L = length(rec_data);
if L == 0
    error('rec_data is empty.');
end

default_nrows = 8;
nrows = min(default_nrows, L);
ncols = floor(L / nrows);
if ncols < 1
    nrows = 1;
    ncols = L;
end

usable = nrows * ncols;
if usable < L
    % 截断多余比特（编码端也会截断）
    rec_data = rec_data(1:usable);
    L = usable;
end

% 反交织
deintlvd_bits = matdeintrlv(rec_data, nrows, ncols);  % 返回按编码端对应顺序的比特

bits_flat = deintlvd_bits(:);  % 列向量

% 把 bits_flat 恢复成编码器输出的比特流（与 f_data_code 配对）
if mod(length(bits_flat), QAMbit) ~= 0
    error('Deinterleaved length not multiple of QAMbit');
end
numSymbols = length(bits_flat) / QAMbit;
temp = reshape(bits_flat, QAMbit, []).';   % numSymbols x QAMbit

coded_stream = temp.';      % QAMbit x numSymbols
coded_stream = coded_stream(:);  % 列向量，应该等于编码端的 coded_bits(1:usable_len)

% ===== Viterbi 译码（硬判决） =====
% 我们不能直接从 trellis 读取 constraintLength 字段（不存在），
% 因此用 numStates 估算约束长度 K (保守估计)：
% 对常见率 1/n: numStates = 2^(K-1) => K = log2(numStates)+1
approx_K = max(1, floor(log2(double(trellis.numStates))) + 1);

% 设置 traceback：通常使用 3~5 倍约束长度；并确保不超过数据长度的一半
traceback = min(5 * approx_K, floor(length(coded_stream) / 2));
if traceback < 1
    traceback = 1;
end

% 执行 Viterbi（hard decision）
y = vitdec(coded_stream, trellis, traceback, 'trunc', 'hard');
y = y(:);
end
