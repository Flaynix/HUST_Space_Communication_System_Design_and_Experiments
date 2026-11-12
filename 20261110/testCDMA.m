%% 演示码分多址原理、远近效应影响及基础的多址干扰消除方法
close all

% A = [1   1    1   1 ];
A = [10 0.1 0.5 1];    % 信号幅度
Tao = [100 200 300 400];  % 信号延迟（码元）


% 生成GPS C/A 码
Code(1,:) = generateCAcode(1);
Code(2,:) = generateCAcode(2);
Code(3,:) = generateCAcode(3);
Code(4,:) = generateCAcode(4);

% 伪码相关性
for m = 1:4
    for n = 1:4
        CCFs{m,n} = ifft(fft(Code(m,:)) .* conj(fft(Code(n,:))));
    end
end
figure;
for m = 1:4
    for n = 1:4
        subplot(4,4,(m-1)*4+n);
        plot(CCFs{m,n});
        xlim([-100 1100]);
        title(['CCFs\{',num2str(m),',',num2str(n),'\}']);
    end
end



% 生成信号
Signal(1,:) = A(1)* circshift(Code(1,:), Tao(1));
Signal(2,:) = A(2)* circshift(Code(2,:), Tao(2));
Signal(3,:) = A(3)* circshift(Code(3,:), Tao(3));
Signal(4,:) = A(4)* circshift(Code(4,:), Tao(4));
SignalRec = sum(Signal) ;
figure;
plot(SignalRec);
xlim([-100 1100]);

% 接收通道相关输出
SigLength = length(Code(1,:));
CorrOut(1,:)  = ifft(fft(SignalRec) .* conj(fft(Code(1,:))))/SigLength;
CorrOut(2,:)  = ifft(fft(SignalRec) .* conj(fft(Code(2,:))))/SigLength;
CorrOut(3,:)  = ifft(fft(SignalRec) .* conj(fft(Code(3,:))))/SigLength;
CorrOut(4,:)  = ifft(fft(SignalRec) .* conj(fft(Code(4,:))))/SigLength;
figure;
for  m =1:4
    subplot(4,1,m);
    plot(CorrOut(m,:));
    xlim([-100 1100]);
    title(['第',num2str(m),'通道的相关输出']);
end

%% 直接估计
for n =1: 4
    [value, index] = max(abs(CorrOut(n,:)));
    AOfEst(n) = value * sign(CorrOut(n,index));
    TaoOfEst(n) = index-1;
end
disp('*********** 直接时延估计 ************');
disp(['设置信号幅度：',num2str(A)]);
disp(['估计信号幅度：',num2str(AOfEst)]);
disp(['设置信号码延迟：',num2str(Tao)]);
disp(['估计信号码延迟：',num2str(TaoOfEst)]);

%% 串行干扰消除
IDestimated     = 1:4 ;   % 已估计的码的序号，按信号功率强度排列
AmpOfestimated  = zeros(1,4);   % 已估计的码的幅度，按信号功率强度排列
TaoOfestimated  = zeros(1,4);   % 已估计的码的相位，按信号功率强度排列
NumOfEstimated  = 0;            % 已估计的码的个数
SignalReduce    = SignalRec;
for m = 1:4
    for n =m: 4
        tempCorr(n,:)  = ifft(fft(SignalReduce) .* conj(fft(Code(IDestimated(n),:))))/SigLength;
        [value, index] = max(abs(tempCorr(n,:)));
        AmpOfestimated(n) = value * sign(tempCorr(n,index));
        TaoOfestimated(n) = index-1;
    end
    [value, index]  = sort(AmpOfestimated(m:4),'descend','ComparisonMethod','abs');
    IDNotestimated  = IDestimated(m:4);
    IDestimated(m:4)     = IDNotestimated(index);
    AmpOfestimated(m:4)  = value;
    TaoOfestimated(m:4)  = TaoOfestimated(m+index-1);
    tempCorr(m,:)        = tempCorr(m+index(1)-1,:);
    SignalReduce =  SignalReduce - AmpOfestimated(m)* circshift(Code(IDestimated(m),:), TaoOfestimated(m));    
    figure;
    subplot(2,1,1);
    plot(tempCorr(m,:));
    xlim([-10,1030]);
    title(['PRN',num2str(IDestimated(m)),'对应信号的相关函数']);
    subplot(2,1,2);
    plot(SignalReduce);
    xlim([-10,1030]);
    title(['消去PRN',num2str(IDestimated(m)),'对应信号之后的残留信号']);
end

figure;
plot(SignalReduce);
title('信号估计残差');

disp('*********** 串行干扰消除 ************');

disp(['设置信号幅度：',num2str(A)]);
AOfEst(IDestimated) = AmpOfestimated;
disp(['估计信号幅度：',num2str(AOfEst)]);
disp(['设置信号码延迟：',num2str(Tao)]);
TaoOfEst(IDestimated) = TaoOfestimated;
disp(['估计信号码延迟：',num2str(TaoOfEst)]);
