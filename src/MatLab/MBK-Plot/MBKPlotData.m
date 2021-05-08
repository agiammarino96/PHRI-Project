%% MBK-Plot

clc
close all
title('MBK-Plot')

xlabel('Virtual Damping B [Nms/rad]')
ylabel('Virtual Stiffness K [Nm/rad]')
grid on
hold on
%Measurement 1
I2 = 0.001;
B2 = [0     0.001   0.005   0.01    0.012   0.015   0.02    0.025   0.027   0.0285];
K2 = [0.00005     0.19    0.39    0.42    0.38    0.32    0.22    0.1     0.05    0];

I3 = 0.01;
B3 = [0     0.001   0.005   0.01    0.05    0.1     0.2     0.25    0.27   ];
K3 = [0.06    0.15    0.95    1.6     4.1     4.1     2.2     0.8     0];


I4 = 0.1;
B4 = [0    0.001    0.05    0.1     0.3     0.5     0.8   0.9   1.1     1.25    1.5     2.0     2.5     2.6 2.7 2.72];
K4 = [0.06    0.1      8.5     16      30      35      41    45    41        39      34     26      10     5   1   0];

plot(B2, K2, '-o', 'LineWidth',1.0)
plot(B3, K3, 'k-o', 'LineWidth',2.0 )
plot(B4, K4, 'k-o', 'LineWidth',2.0 )
%legend('I = 0.001 kgm^2')
%legend('I = 0.001 kgm^2', 'I = 0.01 kgm^2')
legend('I = 0.001 kgm^2', 'I = 0.01 kgm^2', 'I = 0.1 kgm^2')


