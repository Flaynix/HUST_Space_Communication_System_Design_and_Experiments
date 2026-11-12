
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
fs  = 100e6;
% 接收信号的传播延迟 s
TaoOfRec = 100.1 * Tc;
% 接收信号的多普勒频偏 Hz
FdOfRec  = 1000;
% 参考信号的延迟估计 s
TaoOfRef = (0:0.1:200) * Tc;
% 参考信号的多普勒估计 Hz
FdOfRef  = 1000;

% 积分周期 s
Tp  = 1e-3;


%% 仿真相关输出
% 仿真时长
Tsim    = Tp;
t   = 0:1/fs:Tp;
% 接收信号
SigRec  = CAcode(mod(floor((t-TaoOfRec)/Tc),CodeLen)+1) .* exp(1i * 2 *pi* FdOfRec * t);
figure;
plot(t/Tc, real(SigRec),t/Tc, imag(SigRec));
xlabel('时间（Tc）');
title('接收信号波形');
legend('I支路','Q支路');

% 多普勒消除
SigRec1 = SigRec .* exp(-1i * 2 *pi* FdOfRef * t);
figure;
plot(t/Tc, real(SigRec1),t/Tc, imag(SigRec1));
xlabel('时间（Tc）');
title('多普勒消除后的波形');
legend('I支路','Q支路');

% 初始化相关值
ACFs = zeros(size(TaoOfRef));
for m = 1:length(TaoOfRef)
    Tao = TaoOfRef(m);
    RefCode = CAcode(mod(floor((t-Tao)/Tc),CodeLen)+1);
    ACFs(m) = sum(SigRec1 .* RefCode); 
end

% 画出相关输出曲线
DeltaTao = (TaoOfRef - TaoOfRec)/Tc;

figure;
plot(DeltaTao, real(ACFs));
xlabel('时延估计误差（Tc）');
ylabel('I支路相关器输出');

figure;
plot(DeltaTao, imag(ACFs));
xlabel('时延估计误差（Tc）');
ylabel('Q支路相关器输出');

figure;
plot(DeltaTao, abs(ACFs));
xlabel('时延估计误差（Tc）');
ylabel('相关器输出幅值');


