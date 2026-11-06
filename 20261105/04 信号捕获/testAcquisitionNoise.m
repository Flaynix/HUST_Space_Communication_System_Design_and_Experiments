% 仿真在AWGN信道条件下检测概率、虚警概率

% close all;
%% 参数设置
% PRN号
NumOfPRN =10;
% 生成GPS C/A 码
CAcode = generateCAcode(NumOfPRN);
% C/A码码长
CodeLen= 1023;
% 码速率
Rc  = 1.023e6;
% 码元宽度 s
Tc  = 1/Rc;
% 采样率
fs  = Rc;
% 接收信号的传播延迟 s
TaoOfRec = 0 * Tc;
% 接收信号的多普勒频偏 Hz
FdOfRec  = 0;
% 接收信号载噪比
CNRdBs   = 35; 
SNR      = (10.^(CNRdBs/10)/fs);
% 积分周期 s
Tp  = 1e-3;


%% 仿真信号
% 仿真时长
Tsim    = Tp;
t   = 1/fs:1/fs:Tp;
LenSig = length(t);
% 接收信号，加噪声
SigRec  = sqrt(2*SNR)*CAcode(mod(floor((t-TaoOfRec)/Tc),CodeLen)+1) + randn(1,LenSig)+ 1i*randn(1,LenSig);
figure;
plot(t/Tc, real(SigRec),t/Tc, imag(SigRec));
xlabel('时间（Tc）');
title('接收信号波形');
legend('I支路','Q支路');

% 参考信号
SigRef  = CAcode(mod(floor(t/Tc),CodeLen)+1);
figure;
plot(t/Tc, SigRef);
xlabel('时间（Tc）');
title('参考信号波形');
ylim([-1.1 1.1]);
%% 捕获测试
Ntest = 10000;
ACFs = zeros(Ntest,length(SigRec));

for m = 1:Ntest
    SigRec  = sqrt(2*SNR)*CAcode(mod(floor((t-TaoOfRec)/Tc),CodeLen)+1) + randn(1,LenSig)+ 1i*randn(1,LenSig);
    ACFs(m,:) = abs(ifft(fft(SigRec) .* conj(fft(SigRef))));
    waitbar(m/Ntest);
end
figure;
plot(t/Tc,ACFs(m,:))
xlim([-10 1023])
%% 统计检测概率
THs  = 0.1:0.1:5; 
Pd = zeros(size(THs));
Pf = zeros(size(THs));
for m = 1:length(THs)
    Th = THs(m) * sqrt(LenSig);
    % H1
    Pd(m) = nnz(ACFs(:,1)>=Th)/Ntest;
    % H0
    Pf(m) = nnz(ACFs(:,2:LenSig)>=Th)/(Ntest*(LenSig-1));
end

figure;
plot(Pf,Pd)
xlabel('虚警概率');
ylabel('检测概率');
grid on;