
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
TaoOfRec = 100 * Tc;
% 接收信号的多普勒频偏 Hz
FdOfRec  = 1000;
% 积分周期 s
Tp  = 1e-3;


%% 仿真信号
% 仿真时长
Tsim    = Tp;
% 采样时刻
t   = 0:1/fs:Tp;
% 接收信号，无噪声
SigRec  = CAcode(mod(floor((t-TaoOfRec)/Tc),CodeLen)+1) .* exp(1i * 2 *pi* FdOfRec * t);
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
%% 仿真相关输出
DopplerTest = -5000:500:5000; % 多普勒频偏测试点
ACFs = zeros(length(DopplerTest),length(SigRec)); % 相关值矩阵

for m = 1:length(DopplerTest)
    FdOfRef = DopplerTest(m);
    % 多普勒消除
    SigRec1 = SigRec .* exp(-1i * 2 *pi* FdOfRef * t);
    % 基于FFT的码相位并行搜索
    ACFs(m,:) = ifft(fft(SigRec1) .* conj(fft(SigRef)));
end

%%
[xx,yy] = meshgrid(t/Tc,DopplerTest);
figure;
surf(xx,yy,abs(ACFs),'edgecolor','none');
xlabel('码相位(Tc)');
ylabel('多普勒频偏（Hz）');