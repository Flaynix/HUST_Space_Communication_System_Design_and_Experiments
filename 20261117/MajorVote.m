% 多数表决演示
close all;

CAcode1 = generateCAcode(1);
CAcode2 = generateCAcode(2);
CAcode3 = generateCAcode(3);


CAcode4 = generateCAcode(4);
CAcode5 = generateCAcode(5);
CAcode6 = generateCAcode(6);

SI = sign(CAcode1 + CAcode2 + CAcode3);
SQ = sign(CAcode4 + CAcode5 + CAcode6);

S  = sqrt(1/2) * (SI + 1i * SQ);


figure;
subplot(4,1,1)
stairs(CAcode1(1:10),'linewidth',2);
ylabel('分量1');
ylim([-1.2 1.2]);
subplot(4,1,2)
stairs(CAcode2(1:10),'linewidth',2);
ylabel('分量2');
ylim([-1.2 1.2]);
subplot(4,1,3)
stairs(CAcode3(1:10),'linewidth',2);
ylabel('分量3');
ylim([-1.2 1.2]);
subplot(4,1,4)
stairs(SI(1:10),'linewidth',2);
ylabel('合路');
ylim([-1.2 1.2]);

figure;
plot(S,'ro');

Corr1 = mean(S .* CAcode1);
Corr2 = mean(S .* CAcode2);
Corr3 = mean(S .* CAcode3);
Corr4 = mean(S .* CAcode4);
Corr5 = mean(S .* CAcode5);
Corr6 = mean(S .* CAcode6);

A1_6  = abs([Corr1 Corr2 Corr3 Corr4 Corr5 Corr6]);
disp('合路信号中，PRN1-6的幅度分别为：');
disp(num2str(A1_6));
Phi1_6 = angle([Corr1 Corr2 Corr3 Corr4 Corr5 Corr6]) * 180 /pi;
disp('合路信号中，PRN1-6的相位分别为：');
disp(num2str(Phi1_6));
P1_6 = A1_6.^2;
disp('合路信号中，PRN1-6的功率分别为：');
disp(num2str(P1_6));
Eff  = sum(P1_6);
disp('复用效率为：');
disp(num2str(Eff));
