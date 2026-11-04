%% 读取四组环路噪声实验结果
% 先运行 Loop_Bandwidth_Noise_Suppression.slx，再运行本脚本。
% 输入：工作区中的 out；输出：resultsAll、resultsSteady及两张课件形式表。
% 统计结果写入本脚本同目录的result.txt（UTF-8），重复运行覆盖旧结果。
% 本脚本只统计已有数据，不运行或修改Simulink模型，不清空工作区。
% 行标签依据out字段名称，必须确保实际Gain与字段名称标注的带宽一致。
recordSampleTime = 0.01;        % To Workspace的记录间隔，不是噪声发生器采样间隔
recordStartTime = 0;            % 当前模型从0秒开始，To Workspace保留全部数据
steadyStartTime = 50;           % 排除前50秒暂态；建议仿真1000秒获得足够统计样本

assert(exist('out','var')==1, ...
    '工作区没有out。请先运行Loop_Bandwidth_Noise_Suppression.slx。');
assert(isa(out,'Simulink.SimulationOutput') || isstruct(out), ...
    'out应为Simulink.SimulationOutput或包含同名字段的结构体。');
keys = {'1_1Hz','1_01Hz','2_1Hz','2_01Hz'};
LoopOrder = [1;1;2;2];
Bandwidth_Hz = [1;0.1;1;0.1];
% 每行：样本数、起始时间、结束时间、输入噪声功率、误差均方值、RMSE。
allValues = nan(4,6);
steadyValues = nan(4,6);
for k = 1:4
    [noise,tNoise] = readLoggedSignal(out,['Noise_',keys{k}], ...
        recordSampleTime,recordStartTime);
    [err,tErr] = readLoggedSignal(out,['TraErr_',keys{k}], ...
        recordSampleTime,recordStartTime);
    assert(numel(noise)==numel(err) && max(abs(tNoise-tErr))<1e-8, ...
        '第%d组噪声和误差的记录时间不一致，请统一To Workspace设置。',k);
    allValues(k,:) = computeStats(noise,err,tErr);
    keep = tErr >= steadyStartTime;
    if nnz(keep)>=2
        steadyValues(k,:) = computeStats(noise(keep),err(keep),tErr(keep));
    else
        steadyValues(k,1) = 0;
    end
end

%% 输出课件要求的功率与均方根误差
% NoisePower=mean(n.^2)，对应课件rms(Noise).^2。
% ErrorMeanSquare=mean((r-y).^2)，对应课件rms(TraErr).^2。
% RMSE=sqrt(ErrorMeanSquare)；误差含偏置时，均方值不等于去均值后的方差。
columnNames = {'Samples','Start_s','End_s','NoisePower','ErrorMeanSquare','RMSE'};
resultsAll = [table(LoopOrder,Bandwidth_Hz), ...
    array2table(allValues,'VariableNames',columnNames)];
resultsSteady = [table(LoopOrder,Bandwidth_Hz), ...
    array2table(steadyValues,'VariableNames',columnNames)];
fprintf('\n===== 全记录统计（包含启动暂态）=====\n');
disp(resultsAll);
fprintf('NoisePower：输入噪声功率；ErrorMeanSquare：跟踪误差均方值；RMSE：均方根误差。\n');
fprintf('TraErr必须接干净输入减输出，不能接带噪输入减输出。\n');
fprintf('数组时间按记录间隔%.6g秒重建，不使用长度不同的out.tout。\n',recordSampleTime);

fprintf('\n===== 稳态区间统计：t >= %.1f s =====\n',steadyStartTime);
disp(resultsSteady);
% 对应课件第28页：行为带宽1、0.1 Hz，列为一阶、二阶。
errorPowerByBandwidth = table([1;0.1],steadyValues([1,2],5), ...
    steadyValues([3,4],5),'VariableNames', ...
    {'Bandwidth_Hz','FirstOrder_ErrorMeanSquare','SecondOrder_ErrorMeanSquare'});
rmseByBandwidth = table([1;0.1],steadyValues([1,2],6), ...
    steadyValues([3,4],6),'VariableNames', ...
    {'Bandwidth_Hz','FirstOrder_RMSE','SecondOrder_RMSE'});
fprintf('\n课件形式的跟踪误差均方值表：\n');
disp(errorPowerByBandwidth);
fprintf('课件形式的均方根误差表：\n');
disp(rmseByBandwidth);
if any(steadyValues(:,1)==0)
    warning('read_out:TooShort', ...
        ['部分数据未覆盖%.1f秒后的统计区间，稳态结果以NaN表示。' ...
         '请延长模型Stop time（建议1000秒）并重新仿真；不能将全时段值当作稳态值。'], ...
        steadyStartTime);
elseif any(steadyValues(:,3)-steadyValues(:,2)<100)
    warning('read_out:ShortWindow', ...
        '稳态统计窗口不足100秒，尤其0.1 Hz结果的随机偏差可能较大，建议延长仿真。');
end
fprintf('\n参数核对：课件噪声Noise power=1、Sample time=0.01秒。\n');
fprintf('一阶1/0.1 Hz增益：4 / 0.4。\n');
fprintf('二阶1 Hz增益：(3.559986, 2.667925)；0.1 Hz：(0.035600, 0.266792)。\n');
fprintf('仅当实际模型参数正确时，表中的带宽标签及不同组比较才有意义。\n');

%% 保存结果到 result.txt
resultFile = fullfile(fileparts(mfilename('fullpath')),'result.txt');
writeResultText(resultFile,resultsAll,resultsSteady, ...
    errorPowerByBandwidth,rmseByBandwidth,recordSampleTime,steadyStartTime);
fprintf('\n结果已保存至：%s\n',resultFile);

function writeResultText(path,allTable,steadyTable,powerTable,rmseTable,Ts,cutoff)
% 文件句柄仅在函数内部使用，发生异常时也能自动关闭。
[fid,message] = fopen(path,'w','n','UTF-8');
if fid<0
    error('read_out:WriteFailed','无法写入result.txt：%s',message);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'环路带宽对噪声抑制性能的影响\n');
fprintf(fid,'生成时间：%s\n',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf(fid,'模型：Loop_Bandwidth_Noise_Suppression.slx\n');
fprintf(fid,'数组记录间隔：%.9g s；稳态统计起点：%.9g s\n',Ts,cutoff);
fprintf(fid,'NoisePower=mean(Noise.^2)，ErrorMeanSquare=mean(TraErr.^2)，RMSE=sqrt(ErrorMeanSquare)。\n');
fprintf(fid,'TraErr为干净输入减环路输出；行标签依据out字段名，需确认模型实际参数。\n');
fprintf(fid,'\n===== 全记录统计（包含启动暂态）=====\n');
writeNumericTable(fid,allTable);
fprintf(fid,'\n===== 稳态区间统计 =====\n');
writeNumericTable(fid,steadyTable);
fprintf(fid,'\n===== 跟踪误差均方值 =====\n');
writeNumericTable(fid,powerTable);
fprintf(fid,'\n===== 均方根误差 =====\n');
writeNumericTable(fid,rmseTable);
fprintf(fid,'\n字段说明：LoopOrder=环路阶数；Bandwidth_Hz=带宽；Samples=统计样本数；\n');
fprintf(fid,'Start_s/End_s=统计起止时间；NoisePower=输入噪声功率；\n');
fprintf(fid,'ErrorMeanSquare=跟踪误差均方值；RMSE=均方根误差。\n');
if any(steadyTable.Samples==0)
    fprintf(fid,'\n注意：部分记录未覆盖稳态统计区间，对应NaN表示无有效统计结果，不表示零误差。\n');
    fprintf(fid,'请延长仿真时间后重新运行模型和read_out。全记录统计不能代替稳态统计。\n');
elseif any(steadyTable.End_s-steadyTable.Start_s<100)
    fprintf(fid,'\n注意：稳态统计窗口不足100秒，随机偏差可能较大，建议延长仿真。\n');
end
fprintf(fid,'\n课件参数：噪声Noise power=1、Sample time=0.01 s。\n');
fprintf(fid,'一阶1/0.1 Hz增益：4 / 0.4。\n');
fprintf(fid,'二阶1 Hz增益：3.559986、2.667925；0.1 Hz：0.035600、0.266792。\n');
end

function writeNumericTable(fid,T)
% 制表符分隔的纯文本表格，不依赖命令窗口宽度或显示精度。
fprintf(fid,'%s\n',strjoin(T.Properties.VariableNames,sprintf('\t')));
values = table2array(T);
for row=1:size(values,1)
    fprintf(fid,'%.10g',values(row,1));
    fprintf(fid,'\t%.10g',values(row,2:end));
    fprintf(fid,'\n');
end
end

function [x,t] = readLoggedSignal(simOut,name,Ts,t0)
% Array格式没有时间戳，按用户设定的统一记录间隔恢复时间。
try
    value = simOut.(name);
catch
    error('read_out:MissingSignal','out中缺少%s，请检查To Workspace变量名称。',name);
end
if isa(value,'timeseries')
    x = value.Data;
    t = double(value.Time(:));
else
    x = value;
    t = t0+(0:numel(x)-1)'*Ts;
end
assert(isnumeric(x) && isreal(x) && isvector(x) && numel(x)>=2 ...
    && all(isfinite(x(:))),'%s应为至少两个有限实数样本组成的单通道信号。',name);
x = double(x(:));
assert(numel(t)==numel(x) && all(isfinite(t)) && all(diff(t)>0), ...
    '%s的时间戳无效。',name);
assert(max(abs(diff(t)-mean(diff(t))))<1e-8, ...
    '%s不是等间隔采样，请将To Workspace采样时间统一设为0.01秒。',name);
end

function values = computeStats(noise,err,t)
% 等间隔记录的均方统计，不需要额外工具箱。
noisePower = mean(noise.^2);
errorMeanSquare = mean(err.^2);
values = [numel(err),t(1),t(end),noisePower,errorMeanSquare,sqrt(errorMeanSquare)];
end
