clear
clc
% 定义参数
trellis = poly2trellis(7, [171 133]);
M_mod = 4;
data = randi([0 1], 256, 1);

% 编码 + 交织
coded = f_data_code(data, M_mod, trellis);

% 模拟接收端（这里不加噪声）
rec_bits = coded;

% 反交织 + 解码
decoded = f_data_decode(rec_bits, M_mod, trellis);

% 验证误码率
BER = sum(data ~= decoded(1:length(data))) / length(data);
fprintf('BER = %.3e\n', BER);
