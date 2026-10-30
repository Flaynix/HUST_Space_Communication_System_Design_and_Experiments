%% 实验一（1）：多普勒估计正确，1000 Hz
% 输入依次为实验类型、估计频偏、初相误差、输入信噪比、图窗编号。
result1000 = runExperiment('correlation', 1000, 0, Inf, 1);

%% 实验一（2）：多普勒估计为500 Hz
% 500 Hz残余频偏使1 ms相干积分的相关峰降低。
result500 = runExperiment('correlation', 500, 0, Inf, 2);

%% 实验二（1）：初相误差0度，三种DLL鉴别器
resultPhase0 = runExperiment('dll', 1000, 0, Inf, 3);

%% 实验二（2）：初相误差180度，三种DLL鉴别器
% 相干鉴别器反号，包络差及功率差不受整体相位反转影响。
resultPhase180 = runExperiment('dll', 1000, 180, Inf, 4);

%% 实验三（1）：输入信噪比0 dB
resultSNR0 = runExperiment('noise', 1000, 0, 0, 5);

%% 实验三（2）：输入信噪比-10 dB
resultSNRm10 = runExperiment('noise', 1000, 0, -10, 6);

%% 实验三（3）：输入信噪比-20 dB
resultSNRm20 = runExperiment('noise', 1000, 0, -20, 7);

%% 实验三（4）：输入信噪比-50 dB
resultSNRm50 = runExperiment('noise', 1000, 0, -50, 8);

function result = runExperiment(kind, fdEstimate, phaseDegrees, snrDb, figureNumber)
% 核心算法：PRN20生成、接收信号、去载波、码相关、三种鉴别器。
prn = 20;
Rc = 1.023e6;                     % C/A码速率，chip/s
Tc = 1/Rc;
codeLength = 1023;
trueDelay = 30e-3;                % 课件规定的真实传播延迟，s
fdTrue = 1000;
Tp = 1e-3;                       % 课件示例的积分周期，s
spacingChips = 1;                % 超前与滞后总间隔d=1 Tc，各偏移0.5 Tc
fs = 100e6;                      % 对应testACF及testDiscOfDLL示例
if strcmp(kind, 'noise')
    fs = 20e6;                   % 对应课件配套testDiscOfDLLnoise示例
end
N = round(Tp*fs);
t = (0:N-1)/fs;                  % 半开区间[0,Tp)，避免重复计入周期端点
errorChips = (-1000:1000)/10;     % 延迟误差(tauHat-tau)/Tc

% PRN20的G2抽头为4和7；两个寄存器均从全1状态开始。
% G1反馈级：3、10；G2反馈级：2、3、6、8、9、10。
g1 = ones(1,10);
g2 = ones(1,10);
bits = zeros(1,codeLength);
for k = 1:codeLength
    bits(k) = mod(g1(10)+g2(4)+g2(7),2);
    f1 = mod(g1(3)+g1(10),2);
    f2 = mod(sum(g2([2,3,6,8,9,10])),2);
    g1 = [f1,g1(1:9)];
    g2 = [f2,g2(1:9)];
end
assert(bin2dec(sprintf('%d',bits(1:10)))==base2dec('1715',8), ...
    'PRN20生成校验失败。');
caCode = 1-2*bits;
% 与某些生成器相差整体负号不影响同码相关，收发两端采用相同映射。

% r(t)=c(t-tau)*exp(j*(2*pi*fd*t+phi0))+noise。
% 30 ms显式进入码相位；载波传播常相位包含在phi0中。
chipPosition = (t-trueDelay)/Tc;
% 1e-9码片容差只修正浮点运算在整码片边界的舍入误差。
receivedCode = caCode(mod(floor(chipPosition+1e-9),codeLength)+1);
cleanSignal = receivedCode .* exp(1i*(2*pi*fdTrue*t+phaseDegrees*pi/180));
receivedSignal = cleanSignal;
measuredSnrDb = Inf;
if isfinite(snrDb)
    % 与配套示例一致：加入实高斯白噪声，单位信号功率，噪声方差1/SNR。
    % 同一随机序列按不同SNR缩放，使各节独立运行且噪声对比可复现。
    stream = RandStream('mt19937ar','Seed',20261029);
    noise = randn(stream,1,N)/sqrt(10^(snrDb/10));
    receivedSignal = cleanSignal+noise;
    measuredSnrDb = 10*log10(mean(abs(cleanSignal).^2)/mean(abs(noise).^2));
end
carrierRemoved = receivedSignal .* exp(-1i*2*pi*fdEstimate*t);
cleanRemoved = cleanSignal .* exp(-1i*2*pi*fdEstimate*t);

% 扩展扫描范围0.5 Tc，复用相关输出获得E、P、L，避免重复计算。
% ZE(epsilon)=Z(epsilon-0.5)，ZL(epsilon)=Z(epsilon+0.5)。
halfSpacingSteps = round((spacingChips/2)/0.1);
scanChips = (-1000-halfSpacingSteps:1000+halfSpacingSteps)/10;
Z = complex(zeros(size(scanChips)));
Zclean = Z;
for k = 1:numel(scanChips)
    localCode = caCode(mod(floor(chipPosition-scanChips(k)+1e-9),codeLength)+1);
    Z(k) = sum(carrierRemoved .* localCode);  % 与课件相同：累加值，不除以N
    if isfinite(snrDb)
        Zclean(k) = sum(cleanRemoved .* localCode);
    end
end
if ~isfinite(snrDb)
    Zclean = Z;
end
promptIndex = (1:numel(errorChips))+halfSpacingSteps;
P = Z(promptIndex);
E = Z(promptIndex-halfSpacingSteps);
L = Z(promptIndex+halfSpacingSteps);
Eclean = Zclean(promptIndex-halfSpacingSteps);
Lclean = Zclean(promptIndex+halfSpacingSteps);

% 课件第20页三种鉴别器：相干差、非相干包络差、功率差。
% 此处包络差不除以(E+L)，与第12页的归一化包络式有所区别。
Dcoherent = real(E)-real(L);
Denvelope = abs(E)-abs(L);
Dpower = abs(E).^2-abs(L).^2;
Dclean = [real(Eclean)-real(Lclean); abs(Eclean)-abs(Lclean); ...
    abs(Eclean).^2-abs(Lclean).^2];
D = [Dcoherent; Denvelope; Dpower];

result = struct('prn',prn,'sampleRateHz',fs,'sampleCount',N, ...
    'integrationTimeSeconds',Tp,'trueDelaySeconds',trueDelay, ...
    'fdTrueHz',fdTrue,'fdEstimateHz',fdEstimate,'phaseDegrees',phaseDegrees, ...
    'snrDb',snrDb,'measuredSnrDb',measuredSnrDb,'spacingChips',spacingChips, ...
    'delayErrorChips',errorChips,'I',real(P),'Q',imag(P),'magnitude',abs(P), ...
    'normalizedCorrelation',P/N,'early',E,'late',L, ...
    'Dcoherent',Dcoherent,'Denvelope',Denvelope,'Dpower',Dpower, ...
    'cleanDiscriminators',Dclean);
[result.peakMagnitude,peakIndex] = max(abs(P));
result.peakDelayErrorChips = errorChips(peakIndex);
zeroIndex = find(errorChips==0,1);
if ~isfinite(snrDb)
    expected = sum(exp(1i*(2*pi*(fdTrue-fdEstimate)*t+phaseDegrees*pi/180)));
    assert(abs(P(zeroIndex)-expected)<1e-7*N,'对齐点相关校验失败。');
end
fprintf('\n图%d | PRN20 | fs=%.0f MHz | 积分=1 ms | 延迟=30 ms\n',figureNumber,fs/1e6);
fprintf('估计多普勒=%g Hz；初相误差=%g度；SNR=%g dB\n',fdEstimate,phaseDegrees,snrDb);
if strcmp(kind,'correlation')
    fprintf('零误差：I=%.6f，Q=%.6f，幅度=%.6f；归一化幅度=%.6f\n', ...
        real(P(zeroIndex)),imag(P(zeroIndex)),abs(P(zeroIndex)),abs(P(zeroIndex))/N);
else
    fprintf('零误差鉴别输出：相干差=%.6g，包络差=%.6g，功率差=%.6g\n',D(:,zeroIndex));
    if isfinite(snrDb)
        fprintf('本次实测输入SNR=%.3f dB（单次噪声实现）\n',measuredSnrDb);
    end
end

fig = figure(figureNumber);
clf(fig);
set(fig,'Color','w','NumberTitle','off','Position',[60,60,1180,800]);
if strcmp(kind,'correlation')
    values = [real(P);imag(P);abs(P)];
    labels = {'I支路','Q支路','相关幅度'};
    heading = sprintf('实验一 | PRN20 | 多普勒估计 %g Hz | 积分1 ms',fdEstimate);
else
    values = D;
    labels = {'相干超前减滞后','非相干超前减滞后包络','超前功率减滞后功率'};
    if strcmp(kind,'noise')
        heading = sprintf('实验三 | SNR = %g dB | 初相误差0度 | d = 1 T_c',snrDb);
    else
        heading = sprintf('实验二 | 初相误差 %g度 | d = 1 T_c',phaseDegrees);
    end
end
set(fig,'Name',heading);
for row = 1:3
    ax = subplot(3,1,row);
    plot(errorChips,values(row,:),'b-','LineWidth',1.1);
    hold on;
    if isfinite(snrDb)
        plot(errorChips,Dclean(row,:),'k--','LineWidth',1.0);
        if row==1
            legend('含噪声输出','同采样率无噪声基准','Location','best');
        end
    end
    yline(0,':','Color',[0.5,0.5,0.5],'HandleVisibility','off');
    grid on;
    xlabel('延迟估计误差 (\tau_{hat}-\tau)/T_c');
    ylabel('输出（累加值）');
    xlim([-100,100]);
    title([labels{row},'：完整扫描']);
    if strcmp(kind,'correlation')
        ylim([-1.05,1.05]*N);
        if row==3, ylim([-0.05,1.05]*N); end
    elseif ~isfinite(snrDb)
        % 相位0/180度使用相同尺度，直观看出相干曲线反号。
        scale = N;
        if row==3, scale=N^2; end
        ylim([-1.1,1.1]*scale);
    end
    ax.Toolbar.Visible = 'off';
end
sgtitle(heading);
end