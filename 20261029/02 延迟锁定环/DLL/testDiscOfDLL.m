
%% 参数设置
% PRN号
NumOfPRN =1;
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
% 接收信号的初始相位： 弧度
Phi0     = 0 * pi/180;
% 参考信号的延迟估计 s
TaoOfRef = (0:0.1:200) * Tc;
% 参考信号的多普勒估计 Hz
FdOfRef  = 1000;

% 积分周期 s
Tp  = 1e-3;
% 相关器间隔,
d   = 1 * Tc;

%% 仿真相关输出
% 仿真时长
Tsim    = Tp;
t   = 0:1/fs:Tp;
% 接收信号
SigRec  = CAcode(mod(floor((t-TaoOfRec)/Tc),CodeLen)+1) .* exp(1i * 2 *pi* FdOfRec * t) * exp(1i * Phi0);
% 多普勒消除
SigRec1 = SigRec .* exp(-1i * 2 *pi* FdOfRef * t);

% 初始化鉴别器输出
DiscOuts1 = zeros(size(TaoOfRef));  % 相干超前减滞后鉴别器
DiscOuts2 = zeros(size(TaoOfRef));  % 非相干超前减滞后包络鉴别器
DiscOuts3 = zeros(size(TaoOfRef));  % 超前功率减滞后功率鉴别器
for m = 1:length(TaoOfRef)
    Tao = TaoOfRef(m);
    RefCodeEarly = CAcode(mod(floor((t-Tao +d/2)/Tc),CodeLen)+1);
    RefCodeLate  = CAcode(mod(floor((t-Tao -d/2)/Tc),CodeLen)+1);
%     RefCodePrompt= CAcode(mod(floor((t-Tao)/Tc),CodeLen)+1);
    CorrOutEarly = sum(SigRec1 .* RefCodeEarly );
    CorrOutLate  = sum(SigRec1 .* RefCodeLate  );
%     CorrOutPrompt= sum(SigRec1 .* RefCodePrompt);
    DiscOuts1(m)  = (real(CorrOutEarly) - real(CorrOutLate));      % 相干超前减滞后鉴别器
    DiscOuts2(m)  = (abs(CorrOutEarly) - abs(CorrOutLate));        % 非相干超前减滞后包络鉴别器
    DiscOuts3(m)  = (abs(CorrOutEarly).^2 - abs(CorrOutLate).^2);  % 超前功率减滞后功率鉴别器
end

% 画出相关输出曲线
DeltaTao = (TaoOfRef - TaoOfRec)/Tc;

figure;
plot(DeltaTao, real(DiscOuts1));
xlabel('时延估计误差（Tc）');
ylabel('码相位鉴别器输出');
title('相干超前减滞后鉴别器');

figure;
plot(DeltaTao, DiscOuts2);
xlabel('时延估计误差（Tc）');
ylabel('码相位鉴别器输出');
title('超前减滞后包络鉴别器');

figure;
plot(DeltaTao, DiscOuts3);
xlabel('时延估计误差（Tc）');
ylabel('码相位鉴别器输出');
title('超前减滞后功率鉴别器');





