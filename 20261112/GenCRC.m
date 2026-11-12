function CRCout = GenCRC(message, GP, n_CRC)
% message：信息位，MSB 在前；GP：生成多项式非零项的指数列表。
% n_CRC：校验位数；CRCout：由高位到低位排列的校验位行向量。
% 北斗 CRC-24Q 采用零初态、无反射、最终异或为零。

% 将指数列表转换为从高次到低次排列的多项式系数。
message = double(message(:).');
polynomial = zeros(1,n_CRC+1);
polynomial(n_CRC-GP+1) = 1;

% 信息后补 n_CRC 个零，使用异或完成模二长除法。
work = [message, zeros(1,n_CRC)];
for k = 1:numel(message)
    if work(k) == 1
        work(k:k+n_CRC) = xor(work(k:k+n_CRC), polynomial);
    end
end
CRCout = work(end-n_CRC+1:end);
end
