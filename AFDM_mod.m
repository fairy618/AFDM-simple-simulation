% Description: AFDM Modulation
% x: input data vector
% c1, c2: AFDM parameters

function s_blocks = AFDM_mod(X, c1, c2)
% 
% N = size(x,1);
% F = dftmtx(N);
% F = F./norm(F);
% L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
% L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
% A = L2*F*L1;
% out = A'*x;
% 
% end
% 
% 
% function s_blocks = afdm_mod(X, c1, c2)
% 逐列 IDAFT（论文式(1) 的快速实现）
[N, K] = size(X);
s_blocks = zeros(N, K);
for k = 1:K
    s_blocks(:,k) = idaft_col(X(:,k), c1, c2);
end
end

function s = idaft_col(x, c1, c2)
% 单列 IDAFT
[N, K] = size(x); assert(K==1, 'idaft_col 需要 N×1 列向量');
n  = (0:N-1).'; m  = (0:N-1).';
E1 = exp(1j*2*pi*c1*(n.^2));           % 时间域 chirp
E2 = exp(1j*2*pi*c2*(m.^2));           % DAFT 域 chirp
s  = E1 .* (ifft(E2 .* x, [], 1) * sqrt(N));
end