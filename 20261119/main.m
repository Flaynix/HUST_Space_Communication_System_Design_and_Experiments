%% 实验1：BPSK(1) 多径误差包络，BW = 2 MHz
% 直达时延为零；单条反射路径，无噪声。BW 为双边总带宽。
Rc = 1.023e6;
Tc = 1/Rc;
cLight = 299792458;
BW = 2e6;
powerRatioDB = -6;
alpha = 10^(powerRatioDB/20); % 功率比转为幅度比，约 0.5012。
delayChip = (0:0.005:1.3).';
phaseDeg = [0 180];

% 独立等概率 BPSK 矩形码片的理想自相关等于单码片脉冲的自相关。
% 用一个矩形码片求等效相关函数，避免有限伪码互相关引入额外偏差。
samplesPerChip = 1024;
Fs = samplesPerChip*Rc;
Nfft = 2^20;
pulse = zeros(1,Nfft);
pulse(1:samplesPerChip) = 1;
f = [0:Nfft/2-1, -Nfft/2:-1]*Fs/Nfft;
H = double(abs(f) <= BW/2); % 理想矩形低通，通带为 [-BW/2, BW/2]。
X = fft(pulse);
% 频域滤波与相关：IFFT{X*H*conj(X)}，等价于先滤波再相关。
Rperiodic = real(ifft(abs(X).^2 .* H));
lagIndex = -3*samplesPerChip:3*samplesPerChip;
lagChip = lagIndex/samplesPerChip;
R = Rperiodic(mod(lagIndex,Nfft)+1)/Rperiodic(1);
Rinterp = spline(lagChip,R);

% 叠加直达和反射相关函数，寻找相关幅度峰值。
% 测距误差 = 估计距离 - 真实距离；直达时延为零，峰右移对应正误差。
searchChip = -0.6:0.002:0.6;
peakChip = zeros(numel(delayChip),2);
errorMeter = zeros(size(peakChip));
options = optimset('TolX',1e-10,'Display','off');
for k = 1:numel(delayChip)
    for p = 1:2
        echoWeight = alpha*exp(1i*phaseDeg(p)*pi/180);
        corrPower = @(u) abs(ppval(Rinterp,u) + ...
            echoWeight*ppval(Rinterp,u-delayChip(k))).^2;
        [~,iPeak] = max(corrPower(searchChip));
        assert(iPeak>1 && iPeak<numel(searchChip), '峰值落在搜索边界，请扩大搜索范围。');
        % 先网格定位，再在峰值相邻区间细化，避免峰值量化成台阶。
        peakChip(k,p) = fminbnd(@(u) -corrPower(u), ...
            searchChip(iPeak-1),searchChip(iPeak+1),options);
    end
end
errorMeter = peakChip*Tc*cLight;
pathDifference = delayChip*Tc*cLight;

figure(1); clf;
set(gcf,'Name','BPSK(1)：2 MHz','NumberTitle','off','Position',[120 120 850 520]);
plot(pathDifference,errorMeter(:,1),'b-','LineWidth',1.5); hold on;
plot(pathDifference,errorMeter(:,2),'r-','LineWidth',1.5);
yline(0,':','HandleVisibility','off');
grid on; xlim([0 pathDifference(end)]);
xlabel('多径相对路径差 / m'); ylabel('多径测距误差 / m');
title('BPSK(1)，BW = 2 MHz，功率比 = -6 dB');
legend('相移 0°','相移 180°','Location','best');

% 汇总结果，在命令窗口显示并以 UTF-8 合并保存到 result.txt。
resultText = [sprintf('===== BPSK(1)，BW = 2 MHz =====\n'), ...
    sprintf('误差定义：估计距离 - 真实距离（峰右移为正）\n'), ...
    sprintf('功率比 = %.4f；幅度比 = %.6f\n',10^(powerRatioDB/10),alpha), ...
    sprintf('时延范围 = 0~1.3 Tc；路径差范围 = 0~%.3f m\n',pathDifference(end))];
for p = 1:2
    [maxError,iMax] = max(abs(errorMeter(:,p)));
    resultText = [resultText, sprintf( ...
        '相移 %3d°：最大绝对误差 = %.4f m，对应有符号误差 = %+.4f m，相对时延 = %.3f Tc\n', ...
        phaseDeg(p),maxError,errorMeter(iMax,p),delayChip(iMax))]; %#ok<AGROW>
end
fprintf('%s',resultText);
outputDir = fileparts(mfilename('fullpath'));
if isempty(outputDir), outputDir = pwd; end
resultFile = fullfile(outputDir,'result.txt');
% 从已有结果中保留另一节；按带宽替换当前节，防止重复追加。
resultBlocks = {'', ''};
if isfile(resultFile)
    oldText = fileread(resultFile);
    [blockStart,~,bandwidthToken] = regexp(oldText, ...
        '(?m)^===== BPSK\(1\)，BW = (2|24) MHz =====', ...
        'start','end','tokens');
    for b = 1:numel(blockStart)
        if b < numel(blockStart), lastChar = blockStart(b+1)-1;
        else, lastChar = numel(oldText); end
        blockSlot = 1 + strcmp(bandwidthToken{b}{1},'24');
        resultBlocks{blockSlot} = [strtrim(oldText(blockStart(b):lastChar)), newline];
    end
end
resultBlocks{1} = resultText;
combinedText = strjoin(resultBlocks(~cellfun('isempty',resultBlocks)),newline);
[fid,message] = fopen(resultFile,'w','n','UTF-8');
assert(fid ~= -1,'无法创建结果文件：%s',message);
fprintf(fid,'%s',combinedText);
fclose(fid);
fprintf('结果已保存：%s\n',resultFile);

%% 实验2：BPSK(1) 多径误差包络，BW = 24 MHz
% 直达时延为零；单条反射路径，无噪声。BW 为双边总带宽。
Rc = 1.023e6;
Tc = 1/Rc;
cLight = 299792458;
BW = 24e6;
powerRatioDB = -6;
alpha = 10^(powerRatioDB/20); % 功率比转为幅度比，约 0.5012。
delayChip = (0:0.005:1.3).';
phaseDeg = [0 180];

% 独立等概率 BPSK 矩形码片的理想自相关等于单码片脉冲的自相关。
% 用一个矩形码片求等效相关函数，避免有限伪码互相关引入额外偏差。
samplesPerChip = 1024;
Fs = samplesPerChip*Rc;
Nfft = 2^20;
pulse = zeros(1,Nfft);
pulse(1:samplesPerChip) = 1;
f = [0:Nfft/2-1, -Nfft/2:-1]*Fs/Nfft;
H = double(abs(f) <= BW/2); % 理想矩形低通，通带为 [-BW/2, BW/2]。
X = fft(pulse);
% 频域滤波与相关：IFFT{X*H*conj(X)}，等价于先滤波再相关。
Rperiodic = real(ifft(abs(X).^2 .* H));
lagIndex = -3*samplesPerChip:3*samplesPerChip;
lagChip = lagIndex/samplesPerChip;
R = Rperiodic(mod(lagIndex,Nfft)+1)/Rperiodic(1);
Rinterp = spline(lagChip,R);

% 叠加直达和反射相关函数，寻找相关幅度峰值。
% 测距误差 = 估计距离 - 真实距离；直达时延为零，峰右移对应正误差。
searchChip = -0.6:0.002:0.6;
peakChip = zeros(numel(delayChip),2);
errorMeter = zeros(size(peakChip));
options = optimset('TolX',1e-10,'Display','off');
for k = 1:numel(delayChip)
    for p = 1:2
        echoWeight = alpha*exp(1i*phaseDeg(p)*pi/180);
        corrPower = @(u) abs(ppval(Rinterp,u) + ...
            echoWeight*ppval(Rinterp,u-delayChip(k))).^2;
        [~,iPeak] = max(corrPower(searchChip));
        assert(iPeak>1 && iPeak<numel(searchChip), '峰值落在搜索边界，请扩大搜索范围。');
        % 先网格定位，再在峰值相邻区间细化，避免峰值量化成台阶。
        peakChip(k,p) = fminbnd(@(u) -corrPower(u), ...
            searchChip(iPeak-1),searchChip(iPeak+1),options);
    end
end
errorMeter = peakChip*Tc*cLight;
pathDifference = delayChip*Tc*cLight;

figure(2); clf;
set(gcf,'Name','BPSK(1)：24 MHz','NumberTitle','off','Position',[120 120 850 520]);
plot(pathDifference,errorMeter(:,1),'b-','LineWidth',1.5); hold on;
plot(pathDifference,errorMeter(:,2),'r-','LineWidth',1.5);
yline(0,':','HandleVisibility','off');
grid on; xlim([0 pathDifference(end)]);
xlabel('多径相对路径差 / m'); ylabel('多径测距误差 / m');
title('BPSK(1)，BW = 24 MHz，功率比 = -6 dB');
legend('相移 0°','相移 180°','Location','best');

% 汇总结果，在命令窗口显示并以 UTF-8 合并保存到 result.txt。
resultText = [sprintf('===== BPSK(1)，BW = 24 MHz =====\n'), ...
    sprintf('误差定义：估计距离 - 真实距离（峰右移为正）\n'), ...
    sprintf('功率比 = %.4f；幅度比 = %.6f\n',10^(powerRatioDB/10),alpha), ...
    sprintf('时延范围 = 0~1.3 Tc；路径差范围 = 0~%.3f m\n',pathDifference(end))];
for p = 1:2
    [maxError,iMax] = max(abs(errorMeter(:,p)));
    resultText = [resultText, sprintf( ...
        '相移 %3d°：最大绝对误差 = %.4f m，对应有符号误差 = %+.4f m，相对时延 = %.3f Tc\n', ...
        phaseDeg(p),maxError,errorMeter(iMax,p),delayChip(iMax))]; %#ok<AGROW>
end
fprintf('%s',resultText);
outputDir = fileparts(mfilename('fullpath'));
if isempty(outputDir), outputDir = pwd; end
resultFile = fullfile(outputDir,'result.txt');
% 从已有结果中保留另一节；按带宽替换当前节，防止重复追加。
resultBlocks = {'', ''};
if isfile(resultFile)
    oldText = fileread(resultFile);
    [blockStart,~,bandwidthToken] = regexp(oldText, ...
        '(?m)^===== BPSK\(1\)，BW = (2|24) MHz =====', ...
        'start','end','tokens');
    for b = 1:numel(blockStart)
        if b < numel(blockStart), lastChar = blockStart(b+1)-1;
        else, lastChar = numel(oldText); end
        blockSlot = 1 + strcmp(bandwidthToken{b}{1},'24');
        resultBlocks{blockSlot} = [strtrim(oldText(blockStart(b):lastChar)), newline];
    end
end
resultBlocks{2} = resultText;
combinedText = strjoin(resultBlocks(~cellfun('isempty',resultBlocks)),newline);
[fid,message] = fopen(resultFile,'w','n','UTF-8');
assert(fid ~= -1,'无法创建结果文件：%s',message);
fprintf(fid,'%s',combinedText);
fclose(fid);
fprintf('结果已保存：%s\n',resultFile);
