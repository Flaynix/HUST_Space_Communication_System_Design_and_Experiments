%% 实验一：AWGN 中的码相位－多普勒二维相关曲面（C/N0 = 40 dB-Hz）
cfg = experimentConfig();
surfaceResult = acquisitionSurface(cfg);
figSurface = figure('Name', '实验一：带噪声的捕获相关曲面', 'Color', 'w');
plotSurface(axes(figSurface), surfaceResult.codePhaseChips, cfg, ...
    surfaceResult.magnitude, sprintf('PRN %d，C/N_0 = %g dB-Hz，T = %g ms', ...
    cfg.prn, cfg.surfaceCNRdBHz, cfg.surfaceIntegrationTime*1e3));
% 参数框读取配置和实际峰值搜索结果，随参数修改自动更新。
parameterText = sprintf(['设定多普勒频偏：%g Hz\n设定码相位偏移：%g Tc\n' ...
    '估计多普勒频偏：%g Hz\n估计码相位偏移：%.5f Tc'], ...
    cfg.trueDopplerHz, cfg.delayChips, ...
    surfaceResult.estimatedDopplerHz, surfaceResult.estimatedDelayChips);
annotation(figSurface, 'textbox', [0.07 0.76 0.33 0.17], ...
    'String', parameterText, 'Interpreter', 'none', 'FontSize', 11, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'EdgeColor', [0.65 0.65 0.65], 'Margin', 8, 'VerticalAlignment', 'top');
fprintf('真实延迟 %.3f Tc，真实频偏 %.0f Hz。\n', cfg.delayChips, cfg.trueDopplerHz);
fprintf('最大峰：延迟 %.5f Tc，频偏 %.0f Hz。\n', ...
    surfaceResult.estimatedDelayChips, surfaceResult.estimatedDopplerHz);
saveFigure(figSurface, cfg, '01_acquisition_surface');
save(fullfile(cfg.outputDir, '01_acquisition_surface.mat'), 'cfg', 'surfaceResult');

%% 实验二：C/N0 = 30、35、40 dB-Hz 时的 ROC 曲线
cfg = experimentConfig();
probabilityResult = probabilitySimulation(cfg);
rocCNR = [30 35 40];
[~, rocRows] = ismember(rocCNR, probabilityResult.cnrDBHz);
figROC = figure('Name', '实验二：ROC 曲线', 'Color', 'w');
hold on;
rocColors = [0 0.35 0.85; 0 0.6 0.2; 0.9 0.15 0.1];
for k = 1:numel(rocCNR)
    % 主曲线遵循课件：正确码相位统计 Pd，其余码相位统计 Pf。
    row = rocRows(k);
    plot(probabilityResult.pfWrongPhase(row, :), probabilityResult.pd(row, :), ...
        'Color', rocColors(k, :), 'LineWidth', 1.8, ...
        'DisplayName', sprintf('C/N_0 = %d dB-Hz', rocCNR(k)));
end
xlabel('虚警概率 P_f'); ylabel('检测概率 P_d');
title('ROC 曲线（T = 1 ms）');
xlim([0 1]); ylim([0 1]); grid on; legend('Location', 'southeast');
rocTable = table(repelem(rocCNR(:), numel(cfg.thresholdSigma)), ...
    repmat(cfg.thresholdSigma(:), numel(rocCNR), 1), ...
    reshape(probabilityResult.pfWrongPhase(rocRows, :).', [], 1), ...
    reshape(probabilityResult.pd(rocRows, :).', [], 1), ...
    'VariableNames', {'CN0_dBHz', 'Threshold_sigma', 'Pf_wrong_phase', 'Pd'});
saveFigure(figROC, cfg, '02_ROC');
writetable(rocTable, fullfile(cfg.outputDir, '02_ROC.csv'));
save(fullfile(cfg.outputDir, 'probability_data.mat'), 'cfg', 'probabilityResult');
fprintf('\n实验二完成：30、35、40 dB-Hz 三条 ROC 曲线已保存。\n');

%% 实验三：门限为 1σ、2σ、3σ 时，检测概率随 C/N0 的变化
cfg = experimentConfig();
probabilityResult = probabilitySimulation(cfg);
thresholdIndex = arrayfun(@(v) find(abs(cfg.thresholdSigma-v) < 1e-10, 1), [1 2 3]);
pdValues = probabilityResult.pd(:, thresholdIndex);
figPd = figure('Name', '实验三：检测概率与载噪比', 'Color', 'w');
plot(cfg.cnrDBHz, pdValues, '-o', 'LineWidth', 1.6, 'MarkerSize', 4);
xlabel('C/N_0 / (dB-Hz)'); ylabel('检测概率 P_d');
title('检测概率（T = 1 ms）');
legend('1\sigma', '2\sigma', '3\sigma', 'Location', 'southeast');
xlim([25 45]); ylim([0 1]); grid on;
pdTable = table(cfg.cnrDBHz(:), pdValues(:, 1), pdValues(:, 2), pdValues(:, 3), ...
    'VariableNames', {'CN0_dBHz', 'Pd_1sigma', 'Pd_2sigma', 'Pd_3sigma'});
disp(pdTable);
saveFigure(figPd, cfg, '03_Pd_vs_CN0');
writetable(pdTable, fullfile(cfg.outputDir, '03_Pd_vs_CN0.csv'));
save(fullfile(cfg.outputDir, 'probability_data.mat'), 'cfg', 'probabilityResult');

%% 实验四：门限为 1σ、2σ、3σ 时，虚警概率随 C/N0 的变化
cfg = experimentConfig();
probabilityResult = probabilitySimulation(cfg);
thresholdIndex = arrayfun(@(v) find(abs(cfg.thresholdSigma-v) < 1e-10, 1), [1 2 3]);
pfValues = probabilityResult.pfWrongPhase(:, thresholdIndex);
figPf = figure('Name', '实验四：虚警概率与载噪比', 'Color', 'w');
axPf = axes(figPf); hold(axPf, 'on');
pfColors = lines(3);
for k = 1:3
    semilogy(axPf, cfg.cnrDBHz, pfValues(:, k), '-o', ...
        'Color', pfColors(k, :), 'LineWidth', 1.6, 'MarkerSize', 4, ...
        'DisplayName', sprintf('%d\\sigma', k));
end
set(axPf, 'YScale', 'log');
xlabel('C/N_0 / (dB-Hz)'); ylabel('虚警概率 P_f');
title('虚警概率（T = 1 ms）');
xlim([25 45]); ylim([0.005 1]); grid on; legend('Location', 'eastoutside');
pfTable = table(cfg.cnrDBHz(:), pfValues(:, 1), pfValues(:, 2), pfValues(:, 3), ...
    'VariableNames', {'CN0_dBHz', 'Pf_wrong_1sigma', 'Pf_wrong_2sigma', 'Pf_wrong_3sigma'});
disp(pfTable);
saveFigure(figPf, cfg, '04_Pf_vs_CN0');
writetable(pfTable, fullfile(cfg.outputDir, '04_Pf_vs_CN0.csv'));
save(fullfile(cfg.outputDir, 'probability_data.mat'), 'cfg', 'probabilityResult');

% 以下局部函数实现参数配置、FFT 捕获、蒙特卡洛统计和结果保存。
function cfg = experimentConfig()
    cfg.prn = 10;
    cfg.codeRate = 1.023e6;             % C/A 码速率，单位 chip/s
    cfg.surfaceIntegrationTime = 5e-3; % 实验一：5 ms 相干积分，覆盖五个 C/A 码周期
    cfg.surfaceSampleRate = 100e6;     % 与课件无噪声曲面示例一致
    cfg.surfaceCNRdBHz = 40;           % 实验一改为课件要求的 AWGN
    cfg.delayChips = 100;              % 传播延迟 = 100 Tc
    cfg.trueDopplerHz = 1000;
    cfg.dopplerHz = -5000:500:5000;
    cfg.cnrDBHz = 25:1:45;
    cfg.thresholdSigma = (1:50)/10;    % 门限单位：相关输出每个正交支路的噪声标准差
    cfg.trials = 10000;                % 每个载噪比的独立噪声试验次数
    cfg.batchSize = 200;               % 分批计算，避免存储所有试验的完整相关矩阵
    cfg.surfaceSeed = 20261107;
    cfg.probabilitySeed = 20261106;
    cfg.sourceDir = fileparts(mfilename('fullpath'));
    cfg.outputDir = fullfile(cfg.sourceDir, 'results');
    if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
    % 复用老师提供的 C/A 码发生器；核心捕获与概率统计算法均在本文件。
    if ~strcmp(which('generateCAcode'), fullfile(cfg.sourceDir, 'generateCAcode.m'))
        addpath(cfg.sourceDir, '-begin');
    end
end

function result = acquisitionSurface(cfg)
    code = generateCAcode(cfg.prn);
    fs = cfg.surfaceSampleRate;
    n = round(fs * cfg.surfaceIntegrationTime);
    t = (0:n-1)/fs;                    % 500,000 个样点，不重复包含积分终点
    chipPosition = (0:n-1)*cfg.codeRate/fs;
    % 小容差仅修正落在整数码片边界时的浮点舍入。
    localCode = code(mod(floor(chipPosition + 1e-10), 1023)+1);
    receivedCode = code(mod(floor(chipPosition-cfg.delayChips + 1e-10), 1023)+1);
    snrLinear = 10^(cfg.surfaceCNRdBHz/10)/fs;
    stream = RandStream('mt19937ar', 'Seed', cfg.surfaceSeed);
    referenceSpectrum = conj(fft(localCode));
    noise = randn(stream, 1, n) + 1i*randn(stream, 1, n);
    % I、Q 噪声各自方差为 1，总噪声功率为 2；信号功率为 2*SNR。
    received = sqrt(2*snrLinear) * receivedCode .* ...
        exp(1i*2*pi*cfg.trueDopplerHz*t) + noise;
    % 参考码每 1 ms 重复，5 ms 循环相关在码延迟方向也每 1 ms 重复。
    % 所有样点都参与相干积分，仅保留一个码周期的不同延迟候选。
    samplesPerCode = round(fs * numel(code)/cfg.codeRate);
    assert(mod(n, samplesPerCode) == 0, '积分时长必须为完整码周期的整数倍。');
    result.magnitude = zeros(numel(cfg.dopplerHz), samplesPerCode);
    for k = 1:numel(cfg.dopplerHz)
        wiped = received .* exp(-1i*2*pi*cfg.dopplerHz(k)*t);
        correlation = ifft(fft(wiped).*referenceSpectrum);
        % 与示例相同，绘制原始相关累加幅度。
        result.magnitude(k, :) = abs(correlation(1:samplesPerCode));
    end
    [result.peakMagnitude, linearIndex] = max(result.magnitude(:));
    [frequencyIndex, codeIndex] = ind2sub(size(result.magnitude), linearIndex);
    result.codePhaseChips = chipPosition(1:samplesPerCode);
    result.peakCodeIndex = codeIndex;
    result.estimatedDelayChips = chipPosition(codeIndex);
    result.estimatedDopplerHz = cfg.dopplerHz(frequencyIndex);
    result.trueDelayChips = cfg.delayChips;
    result.trueDopplerHz = cfg.trueDopplerHz;
    result.sampleRateHz = fs;
    result.integrationTimeSeconds = cfg.surfaceIntegrationTime;
    result.integrationSampleCount = n;
    result.sigmaRawCorrelation = sqrt(n);
end

function plotSurface(ax, phases, cfg, magnitudeValues, panelTitle)
    [peak, idx] = max(magnitudeValues(:));
    [frequencyIndex, codeIndex] = ind2sub(size(magnitudeValues), idx);
    [~, trueCodeIndex] = min(abs(phases-cfg.delayChips));
    [~, trueFrequencyIndex] = min(abs(cfg.dopplerHz-cfg.trueDopplerHz));
    % 仅显示抽稀，完整数据用于峰值搜索；两个标记所在列额外保留。
    displayIndex = unique([1:10:numel(phases), codeIndex, trueCodeIndex, numel(phases)]);
    surf(ax, phases(displayIndex), cfg.dopplerHz, magnitudeValues(:, displayIndex), 'EdgeColor', 'none');
    hold(ax, 'on');
    peakMarker = plot3(ax, phases(codeIndex), cfg.dopplerHz(frequencyIndex), peak, ...
        'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    trueMarker = plot3(ax, phases(trueCodeIndex), cfg.dopplerHz(trueFrequencyIndex), ...
        magnitudeValues(trueFrequencyIndex, trueCodeIndex), 'kd', 'MarkerSize', 9, 'LineWidth', 1.5);
    legend(ax, [peakMarker trueMarker], {'最大峰', '真实位置'}, 'Location', 'northeast');
    xlabel(ax, '码相位 / T_c'); ylabel(ax, '多普勒频偏 / Hz');
    zlabel(ax, '相关幅度 |R|'); title(ax, panelTitle);
    xlim(ax, [0 1023]); ylim(ax, [cfg.dopplerHz(1), cfg.dopplerHz(end)]);
    view(ax, 40, 35); grid(ax, 'on'); colorbar(ax);
end

function result = probabilitySimulation(cfg)
    % 固定随机种子与完整配置共同作为缓存键，各节均可独立计算。
    persistent cachedConfig cachedResult
    if isequal(cachedConfig, cfg) && ~isempty(cachedResult)
        fprintf('复用相同参数的蒙特卡洛结果。\n');
        result = cachedResult;
        return;
    end
    code = generateCAcode(cfg.prn);
    n = numel(code);                   % fs = Rc，一个码片一个样点，N = 1023
    referenceSpectrum = conj(fft(code));
    codeCorrelation = real(ifft(fft(code).*referenceSpectrum))/sqrt(n);
    amplitude = sqrt(2*10.^(cfg.cnrDBHz/10)/cfg.codeRate);
    nCNR = numel(cfg.cnrDBHz);
    nThreshold = numel(cfg.thresholdSigma);
    detected = zeros(nCNR, nThreshold);
    falseWrongPhase = zeros(nCNR, nThreshold);
    stream = RandStream('mt19937ar', 'Seed', cfg.probabilitySeed);
    fprintf('蒙特卡洛统计：%d 个载噪比，每个 %d 次，T = 1 ms。\n', nCNR, cfg.trials);
    for firstTrial = 1:cfg.batchSize:cfg.trials
        batch = min(cfg.batchSize, cfg.trials-firstTrial+1);
        noise = randn(stream, batch, n) + 1i*randn(stream, batch, n);
        % FFT 沿第 2 维执行，每行对应一次独立试验。
        noiseCorrelation = ifft(fft(noise, [], 2).*referenceSpectrum, [], 2)/sqrt(n);
        for c = 1:nCNR
            % FFT 相关是线性的，先相关噪声再叠加确定信号相关，
            % 与逐次对“信号+噪声”求 FFT 完全等价。
            % 不同 C/N0 共用噪声实现，以减小比较曲线的随机波动。
            correlation = noiseCorrelation + amplitude(c)*codeCorrelation;
            detected(c, :) = detected(c, :) + exceedCounts(abs(correlation(:, 1)), cfg.thresholdSigma);
            falseWrongPhase(c, :) = falseWrongPhase(c, :) + ...
                exceedCounts(abs(correlation(:, 2:end)), cfg.thresholdSigma);
        end
        completed = firstTrial+batch-1;
        if mod(completed, 2000) == 0 || completed == cfg.trials
            fprintf('  已完成 %d / %d 次。\n', completed, cfg.trials);
            drawnow;
        end
    end
    result.cnrDBHz = cfg.cnrDBHz;
    result.thresholdSigma = cfg.thresholdSigma;
    result.pd = detected/cfg.trials;
    result.pfWrongPhase = falseWrongPhase/(cfg.trials*(n-1));
    result.trials = cfg.trials;
    result.samplesPerTrial = n;
    result.sampleRateHz = cfg.codeRate;
    result.sigmaRawCorrelation = sqrt(n);
    result.probabilitySeed = cfg.probabilitySeed;
    result.note = 'Pd: true phase; PfWrongPhase: other phases with signal present. per-cell probabilities.';
    cachedConfig = cfg;
    cachedResult = result;
end

function counts = exceedCounts(values, thresholds)
    % histcounts 的区间左闭右开：反向累加给出 >= 各门限的样本数。
    % 避免为每个门限再次扫描整个相关矩阵。
    bins = histcounts(values(:), [-Inf, thresholds, Inf]);
    counts = fliplr(cumsum(fliplr(bins(2:end))));
end

function saveFigure(fig, cfg, basename)
    set(fig, 'Position', [80 80 1100 680]);
    plotAxes = findall(fig, 'Type', 'axes');
    for k = 1:numel(plotAxes)
        plotAxes(k).Toolbar.Visible = 'off';
    end
    drawnow;
    exportgraphics(fig, fullfile(cfg.outputDir, [basename '.png']), 'Resolution', 180);
    savefig(fig, fullfile(cfg.outputDir, [basename '.fig']));
end
