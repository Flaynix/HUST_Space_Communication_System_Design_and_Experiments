% message 信息位，MSB先发
message = rand(1,576)>0.5 ;
% GP 生成多项式系数
GP = [0,1,3,4,5,6,7,10,11,14,17,18,23,24];
% n_CRC 校验位数
n_CRC =24;
% CRCout 输出的校验位
CRCout = GenCRC(message, GP, n_CRC);   

% ======= 兼容新版MATLAB的改动（替代原来的 crc.generator） =======
h = comm.CRCGenerator('Polynomial', [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1], 'InitialConditions', zeros(1,24), 'FinalXOR', zeros(1,24));
encoded = h(message.'); % 调用生成
% ==================================================================

CRCoutS = encoded(end-n_CRC+1:end).'; % 这行是老师原本代码里提取校验位的

%% 
disp('我的CRC校验位输出：');
disp(num2str(CRCout));
disp('标准CRC校验位输出：');
disp(num2str(CRCoutS));
if CRCout==CRCoutS
    disp('CRC校验位生成正确！');
else
    disp('CRC校验位生成错误！');
end