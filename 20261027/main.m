%% 1. 实验参数与约定
% 伪随机码生成原理与仿真
%
% 统一约定：
% (1) 寄存器数组按课件图中“从左到右”的物理位置编号为1～10。
% (2) 两种结构都向右移位，最右侧寄存器输出；先读取输出，再更新全部寄存器。
% (3) 初态为全1，输出逻辑码 C1、C2；相关运算时映射为 x = 1-2*C。
% (4) 循环相关采用 R12(k)=sum_n x1(n)*x2((n-k) mod N)，不归一化。
%     数学下标 n、k 从0开始，MATLAB数组下标从1开始，R12(k+1)对应延迟k。
%     正k表示将第二个序列循环右移k位后与第一个序列逐点相乘求和。

clearvars;
clc;
n = 10;                              % 寄存器级数
N = 2^n - 1;                         % 理论最大周期：1023码片
initialState = ones(1, n);            % 非零初态：10个1
polyExponents = [10 9 8 6 3 2 0];     % P(X)=X^10+X^9+X^8+X^6+X^3+X^2+1
tapsI = [2 3 6 8 9 10];             % I型：参加反馈异或的寄存器编号

% II型图的多项式标号从左到右递减，与I型图相反。
% 内部项X^9、X^8、X^6、X^3、X^2对应物理寄存器2、3、5、8、9的输入端。
% 对指数e，输入端的数组编号为n-e+1；X^10对应首级直接接回输出。
internalExponents = polyExponents(2:end-1);
tapsII = n - internalExponents + 1;  % 得到[2 3 5 8 9]，不能照抄I型抽头
chipIndex = 0:N-1;
lags = 0:N-1;

fprintf('\n===== 参数设置 =====\n');
fprintf('P(X) = 1 + X^2 + X^3 + X^6 + X^8 + X^9 + X^10\n');
fprintf('寄存器级数 n = %d，理论周期 N = %d 码片，初始状态为全1。\n', n, N);

%% 2. 通过I型结构生成C1
% I型为外部反馈结构：多个抽头先异或，再送入第1级。
% feedback = r(2) XOR r(3) XOR r(6) XOR r(8) XOR r(9) XOR r(10)
% new_r = [feedback, old_r(1), ..., old_r(9)]。
% mod(sum(...),2)就是多个二进制位的异或，即模二和。

assert(exist('tapsI', 'var') == 1, '缺少寄存器参数。');
regI = initialState;
C1 = zeros(1, N);
stateHistoryI = zeros(N, n);          % 记录输出前的状态，用于验证实际周期
for chip = 1:N
    stateHistoryI(chip, :) = regI;
    C1(chip) = regI(end);            % 读取更新前的最右侧寄存器
    feedbackI = mod(sum(regI(tapsI)), 2);
    regI = [feedbackI, regI(1:end-1)]; % 同一时钟沿并行更新，避免串行覆盖旧状态
end

% N个状态两两不同，且第N次更新后回到初态，因此实际周期恰好为N。
periodI = size(unique(stateHistoryI, 'rows'), 1);
assert(periodI == N && isequal(regI, initialState), 'I型寄存器未达到1023位周期。');
assert(sum(C1) == 2^(n-1), 'C1未满足m序列的平衡性。');
x1 = 1 - 2*C1;                      % 逻辑电平0 -> +1，逻辑电平1 -> -1

fprintf('\n===== I型结构：C1 =====\n');
fprintf('实测周期 = %d；1的个数 = %d；0的个数 = %d。\n', periodI, sum(C1), N-sum(C1));
fprintf('前32位：%s\n', sprintf('%d', C1(1:32)));

figure(1); clf;
set(gcf, 'Name', '实验节点2：I型结构生成C1', 'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1100 700]);
stairs(0:N, [x1 x1(end)], 'LineWidth', 0.8);
grid on; xlim([0 N]); ylim([-1.5 1.5]); yticks([-1 1]);
xlabel('码片序号 n'); ylabel('双极性码 x_1'); title('C_1波形：一个完整周期，0映射为+1、1映射为-1');
drawnow;

%% 3. 通过II型结构生成C2
% II型为内部反馈结构：最后一级输出接回首级，同时送到各内部异或点。
% 每次先保存整个旧状态，再构造新状态，所有异或都使用同一个旧输出。
% 设旧状态为[r1,...,r10]，则课件第16页II型图对应：
% new_r = [r10, r1 XOR r10, r2 XOR r10, r3, r4 XOR r10, ...
%          r5, r6, r7 XOR r10, r8 XOR r10, r9]。
% 注意：II型与I型的寄存器坐标含义不同，全1初态不代表输出码相位相同。

assert(exist('tapsII', 'var') == 1, '缺少寄存器参数。');
regII = initialState;
C2 = zeros(1, N);
stateHistoryII = zeros(N, n);
for chip = 1:N
    stateHistoryII(chip, :) = regII;
    oldStateII = regII;
    outputII = oldStateII(end);
    C2(chip) = outputII;
    regII = [outputII, oldStateII(1:end-1)]; % 先构造直接回送与普通移位
    regII(tapsII) = mod(regII(tapsII) + outputII, 2); % 再加入内部反馈异或
end

periodII = size(unique(stateHistoryII, 'rows'), 1);
assert(periodII == N && isequal(regII, initialState), 'II型寄存器未达到1023位周期。');
assert(sum(C2) == 2^(n-1), 'C2未满足m序列的平衡性。');
x2 = 1 - 2*C2;

fprintf('\n===== II型结构：C2 =====\n');
fprintf('实测周期 = %d；1的个数 = %d；0的个数 = %d。\n', periodII, sum(C2), N-sum(C2));
fprintf('前32位：%s\n', sprintf('%d', C2(1:32)));

figure(2); clf;
set(gcf, 'Name', '实验节点3：II型结构生成C2', 'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1100 700]);
stairs(0:N, [x2 x2(end)], 'LineWidth', 0.8);
grid on; xlim([0 N]); ylim([-1.5 1.5]); yticks([-1 1]);
xlabel('码片序号 n'); ylabel('双极性码 x_2'); title('C_2波形：一个完整周期，0映射为+1、1映射为-1');
drawnow;

%% 4. 计算并绘制C1、C2的循环自相关
% 循环自相关：Rii(k) = sum_{n=0}^{N-1} xi(n)*xi((n-k) mod N)。
% 按照定义依次执行循环移位、逐点相乘和求和。
% mod保证末尾与开头相接，没有补零，因此得到循环相关。
% 不使用线性相关函数xcorr，也不直接对0/1逻辑码做相关。
% 未归一化理论值：Rii(0)=N，Rii(k)=-1（k=1,...,N-1）。
% rho11、rho22为归一化自相关结果；理论值为1和-1/N。

assert(exist('C1', 'var') == 1 && exist('C2', 'var') == 1, ...
    '缺少序列C1或C2。');
x1 = 1 - 2*C1;
x2 = 1 - 2*C2;
R11 = zeros(1, N);
R22 = zeros(1, N);
for k = lags
    shiftedIndex = mod(chipIndex - k, N) + 1; % 与循环右移k位等价
    R11(k+1) = sum(x1 .* x1(shiftedIndex)); % C1 循环自相关
    R22(k+1) = sum(x2 .* x2(shiftedIndex)); % C2 循环自相关
end
rho11 = R11 / N;
rho22 = R22 / N;
assert(R11(1) == N && all(R11(2:end) == -1), 'C1自相关不符合m序列理论。');
assert(R22(1) == N && all(R22(2:end) == -1), 'C2自相关不符合m序列理论。');

fprintf('\n===== 循环自相关（未归一化）=====\n');
fprintf('C1：R11(0) = %d，其余延迟的相关值全部为-1。\n', R11(1));
fprintf('C2：R22(0) = %d，其余延迟的相关值全部为-1。\n', R22(1));
fprintf('归一化后：主峰为1，旁瓣为-1/N = %.9f。\n', -1/N);

figure(3); clf;
set(gcf, 'Name', '实验节点4：循环自相关', 'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1100 700]);
subplot(2, 1, 1);
plot(lags, R11, 'b-', 'LineWidth', 1); hold on;
plot(0, R11(1), 'ro', 'MarkerFaceColor', 'r');
grid on; xlim([-10 N]); ylim([-50 N+100]);
xlabel('循环延迟 k / 码片'); ylabel('R_{11}(k)'); title('C_1循环自相关：主峰1023');
subplot(2, 1, 2);
plot(lags, R22, 'b-', 'LineWidth', 1); hold on;
plot(0, R22(1), 'ro', 'MarkerFaceColor', 'r');
grid on; xlim([-10 N]); ylim([-50 N+100]);
xlabel('循环延迟 k / 码片'); ylabel('R_{22}(k)'); title('C_2循环自相关：主峰1023');
drawnow;

%% 5. 计算并绘制C1与C2的循环互相关
% R12(k) = sum_n x1(n)*x2((n-k) mod N)。
% 同时计算反向互相关R21，验证R21(k)=R12((-k) mod N)。
% 如果两序列仅差一个循环移位，互相关将在对齐处为N，其余处为-1。
% 横轴显式使用lags=0:N-1，不能把MATLAB的峰值数组下标直接当成相位差。

assert(exist('C1', 'var') == 1 && exist('C2', 'var') == 1, ...
    '缺少序列C1或C2。');
x1 = 1 - 2*C1;
x2 = 1 - 2*C2;
R12 = zeros(1, N);
R21 = zeros(1, N);
for k = lags
    shiftedIndex = mod(chipIndex - k, N) + 1;
    R12(k+1) = sum(x1 .* x2(shiftedIndex)); % C1 与 C2 循环互相关
    R21(k+1) = sum(x2 .* x1(shiftedIndex)); % 反向循环互相关
end
[peak12, peakIndex12] = max(R12);
[peak21, peakIndex21] = max(R21);
delay12 = lags(peakIndex12);          % MATLAB数组下标减1，得到实际循环延迟
delay21 = lags(peakIndex21);
sideMask12 = lags ~= delay12;
assert(peak12 == N && all(R12(sideMask12) == -1), '互相关未满足单峰两值特性。');
assert(isequal(R21, R12(mod(-lags, N)+1)), '正反向互相关未满足循环对称关系。');
rho12 = R12 / N;

fprintf('\n===== 循环互相关（未归一化）=====\n');
fprintf('R12峰值 = %d，峰值对应延迟k = %d码片（数组下标为%d）。\n', ...
    peak12, delay12, peakIndex12);
fprintf('R21峰值 = %d，峰值对应延迟k = %d码片。\n', peak21, delay21);
fprintf('除唯一主峰外，其余互相关值全部为-1。\n');

figure(4); clf;
set(gcf, 'Name', '实验节点5：循环互相关', 'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1100 700]);
subplot(2, 1, 1);
plot(lags, R12, 'b-', 'LineWidth', 1.2); hold on;
plot(delay12, peak12, 'ro', 'MarkerFaceColor', 'r');
text(delay12+15, peak12, sprintf('k=%d，R_{12}=%d', delay12, peak12));
grid on; xlim([0 N-1]); ylim([-50 N+150]);
xlabel('C_2循环右移量 k / 码片'); ylabel('R_{12}(k)'); title('C_1与C_2的循环互相关');
subplot(2, 1, 2);
nearPeak = max(0, delay12-12):min(N-1, delay12+12);
% 横轴仅显示主峰左右各12个码片，纵轴与全周期图一致，完整保留主峰。
stem(nearPeak, R12(nearPeak+1), 'filled', 'MarkerSize', 3); hold on;
plot(delay12, peak12, 'ro', 'MarkerFaceColor', 'r');
text(delay12+0.5, peak12, sprintf('k=%d，R_{12}=%d', delay12, peak12));
grid on; xlim([nearPeak(1)-1 nearPeak(end)+1]); ylim([-50 N+150]);
xlabel('C_2循环右移量 k / 码片'); ylabel('R_{12}(k)');
title('循环互相关：峰值附近的横轴放大图');
drawnow;

%% 6. 相位差分析与逐码片对齐验证
% 本节用实际移位及逐位比较验证结论，避免只根据曲线观察相位差。
% 按本程序的定义，应得到delay12=320、delay21=703，两者之和为1023。
% C2比C1超前320码片：C2(n)=C1((n+320) mod 1023)。
% 所以将C2循环右移320位后得到C1；等价地，循环左移703位后得到C1。
% 320是“码相位差”，单位为码片。课件未给码片周期Tc，故不换算为秒。
% 若交换互相关的序列顺序，或改用(n+k)定义，峰值将位于703而非320。

assert(exist('R12', 'var') == 1 && exist('C1', 'var') == 1 && exist('C2', 'var') == 1, ...
    '缺少序列C1、C2或循环互相关结果R12。');
[peakValue, peakIndex] = max(R12);
phaseAdvance = peakIndex - 1;
equivalentLeftShift = mod(N-phaseAdvance, N);
C2Aligned = circshift(C2, [0 phaseAdvance]); % 正数表示对行向量循环右移
mismatchBefore = sum(C1 ~= C2);
mismatchAfter = sum(C1 ~= C2Aligned);
assert(peakValue == N && mismatchAfter == 0, '按互相关峰值移位后仍存在不一致码片。');
assert(isequal(C2, circshift(C1, [0 -phaseAdvance])), '相位超前关系验证失败。');

fprintf('\n===== 相位差及对齐验证 =====\n');
fprintf('C2相对C1超前 %d 码片（按一个周期取模）。\n', phaseAdvance);
fprintf('C2循环右移 %d 位，或循环左移 %d 位，均与C1完全重合。\n', ...
    phaseAdvance, equivalentLeftShift);
fprintf('对齐前不一致码片数 = %d/%d。\n', mismatchBefore, N);
fprintf('对齐后不一致码片数 = %d/%d。\n', mismatchAfter, N);
fprintf('结论：两种结构产生同一m序列的不同循环相位；全1寄存器初态不保证码相位相同。\n');
fprintf('全部实验计算完成。\n');
