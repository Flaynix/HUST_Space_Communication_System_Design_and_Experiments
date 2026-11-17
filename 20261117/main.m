%% 实验一：六分量多数表决恒包络复用
% 六路独立等概率二进制序列，映射为 0 -> +1、1 -> -1。
rng(20261117, 'twister');
L = 1e6;
bits = randi([0 1], 6, L);
code = 1 - 2*bits;

% 前三路位于 I 支路，后三路位于 Q 支路；除以 sqrt(2) 使功率为 1。
SI = sign(sum(code(1:3,:), 1));
SQ = sign(sum(code(4:6,:), 1));
S = (SI + 1i*SQ)/sqrt(2);
powerMV = mean(abs(S).^2);

% 三路多数表决的代数式为 (c1+c2+c3-c1*c2*c3)/2。
% 每路有用分量幅度为 1/(2*sqrt(2))，六路总效率为 6/8 = 75%。
corrTheoryMV = [ones(3,1); 1i*ones(3,1)]/(2*sqrt(2));
corrSimMV = code*S.'/L;  % 非共轭转置，保留 Q 支路的正相位。
ampTheoryMV = abs(corrTheoryMV);
ampSimMV = abs(corrSimMV);
phaseTheoryMV = angle(corrTheoryMV)*180/pi;
phaseSimMV = angle(corrSimMV)*180/pi;
effTheoryMV = sum(ampTheoryMV.^2);
effSimMV = sum(ampSimMV.^2)/powerMV;

resultMV = table((1:6).', ampTheoryMV, ampSimMV, ...
    phaseTheoryMV, phaseSimMV, ...
    'VariableNames', {'Component','AmplitudeTheory','AmplitudeSim', ...
    'PhaseTheory_deg','PhaseSim_deg'});
% 同一份结果文本同时显示在命令窗口，并以 UTF-8 覆盖保存。
resultText = [sprintf('===== 多数表决：理论与仿真 =====\n'), ...
    sprintf('样本数 = %d；随机种子 = 20261117\n', L), ...
    evalc('disp(resultMV);'), ...
    sprintf('平均功率 = %.6f；最大包络偏差 = %.3g\n', powerMV, max(abs(abs(S)-1))), ...
    sprintf('理论效率 = %.4f%%；仿真效率 = %.4f%%\n', 100*effTheoryMV, 100*effSimMV)];
fprintf('%s', resultText);
outputDir = fileparts(mfilename('fullpath'));
if isempty(outputDir), outputDir = pwd; end
resultFile = fullfile(outputDir, 'result_majority.txt');
[fid, message] = fopen(resultFile, 'w', 'n', 'UTF-8');
assert(fid ~= -1, '无法创建结果文件：%s', message);
fprintf(fid, '%s', resultText);
fclose(fid);
fprintf('结果已保存：%s\n', resultFile);

figure(1); clf;
set(gcf, 'Name', '六分量多数表决', 'NumberTitle', 'off', 'Position', [100 100 1000 680]);
tiledlayout(2,2, 'TileSpacing','compact');
nexttile;
plot(real(S(1:2000)), imag(S(1:2000)), '.', 'MarkerSize',12);
axis equal; xlim([-1.2 1.2]); ylim([-1.2 1.2]); grid on;
xlabel('I'); ylabel('Q'); title('星座图');
nexttile;
plot(1:100, abs(S(1:100)), 'LineWidth',1.5);
ylim([0 1.2]); grid on; xlabel('码片序号'); ylabel('|S|'); title('合路信号包络');
nexttile;
bar(1:6, [ampTheoryMV ampSimMV]); grid on;
xlabel('分量编号'); ylabel('幅度'); title('分量幅度'); legend('理论','仿真', 'Location','best');
nexttile;
bar(1:6, [phaseTheoryMV phaseSimMV]); grid on;
xlabel('分量编号'); ylabel('相位 / °'); title('分量相位'); legend('理论','仿真', 'Location','best');

%% 实验二：基于老师查找表的六分量 POCET 恒包络复用
% POCET1 提供 64 种输入组合及其相位（单位为度），捕获其打印文本。
teacherTableText = evalc('POCET1;');
phaseLUT = exp(1i*ListPhi*pi/180);

% 理论分析：六路独立等概率输入的 64 种组合等概率出现。
% ListCode 的列顺序与二进制从最高位到最低位的顺序一致。
corrTheoryPOCET = ListCode.'*phaseLUT/64;
ampTheoryPOCET = abs(corrTheoryPOCET);
phaseTheoryPOCET = angle(corrTheoryPOCET)*180/pi;
effTheoryPOCET = sum(ampTheoryPOCET.^2);

% 与实验一采用相同种子和样本数，保持随机输入一致；本节可独立计算。
rng(20261117, 'twister');
L = 1e6;
bits = randi([0 1], 6, L);
code = 1 - 2*bits;
% 第一分量为最高位；十进制组合值加 1 得到 MATLAB 的查表行号。
rowIndex = (2.^(5:-1:0))*bits + 1;
S = reshape(phaseLUT(rowIndex), 1, []);
powerPOCET = mean(abs(S).^2);
corrSimPOCET = code*S.'/L;
ampSimPOCET = abs(corrSimPOCET);
phaseSimPOCET = angle(corrSimPOCET)*180/pi;
effSimPOCET = sum(ampSimPOCET.^2)/powerPOCET;

resultPOCET = table((1:6).', ampTheoryPOCET, ampSimPOCET, ...
    phaseTheoryPOCET, phaseSimPOCET, ...
    'VariableNames', {'Component','AmplitudeTheory','AmplitudeSim', ...
    'PhaseTheory_deg','PhaseSim_deg'});
% 两个实验分别保存，独立运行其中一节不会覆盖另一节的结果。
resultText = [sprintf('===== POCET：理论与仿真 =====\n'), ...
    sprintf('样本数 = %d；随机种子 = 20261117\n', L), ...
    evalc('disp(resultPOCET);'), ...
    sprintf('平均功率 = %.6f；最大包络偏差 = %.3g\n', powerPOCET, max(abs(abs(S)-1))), ...
    sprintf('理论效率 = %.4f%%；仿真效率 = %.4f%%\n', 100*effTheoryPOCET, 100*effSimPOCET)];
fprintf('%s', resultText);
outputDir = fileparts(mfilename('fullpath'));
if isempty(outputDir), outputDir = pwd; end
resultFile = fullfile(outputDir, 'result_pocet.txt');
[fid, message] = fopen(resultFile, 'w', 'n', 'UTF-8');
assert(fid ~= -1, '无法创建结果文件：%s', message);
fprintf(fid, '%s', resultText);
fclose(fid);
fprintf('结果已保存：%s\n', resultFile);

figure(2); clf;
set(gcf, 'Name', '六分量 POCET', 'NumberTitle', 'off', 'Position', [140 100 1000 680]);
tiledlayout(2,2, 'TileSpacing','compact');
nexttile;
plot(real(S(1:2000)), imag(S(1:2000)), '.', 'MarkerSize',12);
axis equal; xlim([-1.2 1.2]); ylim([-1.2 1.2]); grid on;
xlabel('I'); ylabel('Q'); title('星座图');
nexttile;
plot(1:100, abs(S(1:100)), 'LineWidth',1.5);
ylim([0 1.2]); grid on; xlabel('码片序号'); ylabel('|S|'); title('合路信号包络');
nexttile;
bar(1:6, [ampTheoryPOCET ampSimPOCET]); grid on;
xlabel('分量编号'); ylabel('幅度'); title('分量幅度'); legend('理论','仿真', 'Location','best');
nexttile;
bar(1:6, [phaseTheoryPOCET phaseSimPOCET]); grid on;
xlabel('分量编号'); ylabel('相位 / °'); title('分量相位'); legend('理论','仿真', 'Location','best');
