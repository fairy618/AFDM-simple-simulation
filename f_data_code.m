% 我需要实现一个matlab的子函数：function y = f_data_code(data, M_mod, trellis)
% 其中data是需要编码的二进制数据（列向量）；
% M_mod是QAM调制阶数，例如4, 16, 64等（暂时不考虑非标准的，认为我输入的就是2、4、16、64等等）
% trellis是我生成的卷积码的trellis结构。
% 我需要你做的就是，补充函数，包括：
% 1、将数据卷积编码
% 2、将编码的数据重新整形
% 3、将整形后的数据进行交织

function y = f_data_code(data, M_mod, trellis)
% f_data_code : 卷积编码 + 整形 + 交织（按块交织，鲁棒）
% 输入：
%   data    : 二进制列向量
%   M_mod   : QAM阶数 (e.g. 4,16,64)
%   trellis : 卷积码trellis
% 输出：
%   y : 交织后的比特列向量

% 1) 卷积编码
coded_bits = convenc(data(:), trellis);   % 列向量

% 2) QAMbit（用于整形检查）
QAMbit = log2(M_mod);
if QAMbit ~= round(QAMbit)
    error('M_mod must be power of 2');
end

% 3) 整形为按符号组织（每列一个符号的比特）
% 如果长度不是 QAMbit 的整数倍，截断多余比特
len = length(coded_bits);
numSymbols = floor(len / QAMbit);
if numSymbols == 0
    error('Coded bit length too short for the given M_mod.');
end
usable_len = numSymbols * QAMbit;
if usable_len < len
    coded_bits = coded_bits(1:usable_len);
end
bit_matrix = reshape(coded_bits, QAMbit, []).';  % numSymbols x QAMbit

% 4) 我们做按比特流的块交织（block interleaver）
% 将 bit_matrix 展平为按符号顺序的比特向量：
bits_flat = bit_matrix.';    % QAMbit x numSymbols
bits_flat = bits_flat(:);    % 列向量，按每符号 QAMbit 比特顺序排列

L = length(bits_flat);       % 总比特数

% 交织参数自动选择（保证不会出现 ncols = 0）
default_nrows = 8;
nrows = min(default_nrows, L);    % 行数不超过总长度
ncols = floor(L / nrows);
if ncols < 1
    nrows = 1;
    ncols = L;
end
usable = nrows * ncols;
if usable < L
    % 截断末尾多余的比特以配合整块交织
    bits_flat = bits_flat(1:usable);
    L = usable;
end

% 对每个 block（这里 block 就是整段 bits_flat，因为我们确保 L = nrows*ncols）
% 如果总长度很大，也可以分多块，这里直接一次性交织整个 bits_flat
intlvd_bits = matintrlv(bits_flat, nrows, ncols);

% 输出保持为列向量（交织后的比特流）
y = intlvd_bits(:);

end
