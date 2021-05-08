close all
clear all
clc

%% Calibration Force Sensor NEGATIVE
F=-[0 0.164 0.197 0.511 0.737]*9.81;
noise_F=0.01;
V_out_F_ww=[0.08 0.54 0.62 1.42 2.09];

figure(1)
plot(F,V_out_F_ww,'o')
hold on
xlabel('Force [N]')
ylabel('Voltage [V]')
title('Calibration Force Sensor')

n=length(F);
x=F;
y=V_out_F_ww;
x_mean=mean(x);
y_mean=mean(y);
num_beta=0;
den_beta=0;
for i=1:n
    num_beta=num_beta+(x(i)-x_mean)*(y(i)-y_mean);
    den_beta=den_beta+(x(i)-x_mean)^2;
end
beta=num_beta/den_beta;
alpha=y_mean-beta*x_mean;
y=beta*x+alpha;
plot(x,y,'-r','LineWidth',1)
%  x=(y-alpha)/beta
%  alpha=0.0838
%  beta=0.2744
%  Force is positive in -x
%% Calibration Force Sensor POSITIVE

F=[0 0.164 0.197 0.511 0.737]*9.81;
V_out_F_ww=[0.08 -0.41 -0.56 -1.31 -1.94];


plot(F,V_out_F_ww,'o')



n=length(F);
x=F;
y=V_out_F_ww;
x_mean=mean(x);
y_mean=mean(y);
num_beta=0;
den_beta=0;
for i=1:n
    num_beta=num_beta+(x(i)-x_mean)*(y(i)-y_mean);
    den_beta=den_beta+(x(i)-x_mean)^2;
end
beta_1=num_beta/den_beta;
alpha_1=y_mean-beta*x_mean;
y=beta_1*x+alpha_1;
plot(x,y,'-k','LineWidth',1)
legend('Measures Neg','Linear Fit Neg','Meas Pos','Linear Fit Pos')

%% Params useful
r_p1	= 0.075; % [m]          CAD
r_m     = 0.005; % [m]          CAD
r_p2    = 0.065; % [m]          CAD
Jac=r_m/r_p1*r_p2; % [m]
rp1_rm=r_p1/r_m;

%% Transf function DeltaV-->Vel

V=[120 34 26 20 15 9 2.3 -3.5 -9.8 -16 -22.5 -29 -38 -45 -51 -59 -64 -72 -132]*(-rp1_rm); % [deg/s]
DeltaV=[-1 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.10 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 1];      % [V]  

figure(2)
plot(DeltaV,V,'o')
hold on

n=length(DeltaV);
x=DeltaV;
y=V;
x_mean=mean(x);
y_mean=mean(y);
num_beta=0;
den_beta=0;
for i=1:n
    num_beta=num_beta+(x(i)-x_mean)*(y(i)-y_mean);
    den_beta=den_beta+(x(i)-x_mean)^2;
end
beta_2=num_beta/den_beta;
alpha_2=y_mean-beta_2*x_mean;
y=beta_2*x+alpha_2; % Velocity=beta*Voltage+alpha
plot(x,y,'-k','LineWidth',1)
title('Transfer function amplifier')
xlabel('\Delta V [V]')
ylabel('Motor Velocity [°/s]')
legend('Measure','Fitting')
% Voltage=(Vel-alpha)/beta



%% Saturation Force Vel Acc for experiment transparency

Asat=readtable('Accsat.xlsx');
Fsat=readtable('Fsat.xlsx');
Vsat=readtable('Vsat.xlsx');


A_sat=table2array(Asat)*pi/180;
F_sat=table2array(Fsat);
V_sat=table2array(Vsat)*pi/180;

Asaturation=max(abs(A_sat));
Fsaturation=max(abs(F_sat));
Vsaturation=max(abs(V_sat));
Tsaturation=Fsaturation*r_p2;

% Fsaturation=1N
% Asaturation=539.9220 rad/s^2
% Vsaturation= 20.0651 rad/s



%% First test (transparency planes)
Ameas=readtable('AccMeas.xlsx');
Fmeas=readtable('Fmeas.xlsx');
Vmeas=readtable('Vmeas.xlsx');

Ameas_adm=readtable('AccMeasadm.xlsx');
Fmeas_adm=readtable('Fmeasadm.xlsx');
Vmeas_adm=readtable('Vmeasadm.xlsx');


A_adm=table2array(Ameas_adm)*pi/180;
F_adm=table2array(Fmeas_adm);
V_adm=table2array(Vmeas_adm)*pi/180;

Data_adm=[A_adm V_adm F_adm*r_p2];
N=length(A_adm);
B_adm=[ones(N,1) Data_adm(:,1:2)]\Data_adm(:,3);
a_adm=B_adm(2);
b_adm=B_adm(3);
c_adm=B_adm(1);
x_adm=[];
y_adm=[];
z_adm=[];
for i=-300:1:300
    for j=-10:1/30:10
        x_adm(end+1)=i;
        y_adm(end+1)=j;
        z_adm(end+1)=a_adm*i+b_adm*j+c_adm;
    end
end

Ameas_Imp=readtable('AccImp.xlsx');
Fmeas_Imp=readtable('FImp.xlsx');
Vmeas_Imp=readtable('VImp.xlsx');


A_Imp=table2array(Ameas_Imp)*pi/180;
F_Imp=table2array(Fmeas_Imp);
V_Imp=table2array(Vmeas_Imp)*pi/180;

Data_Imp=[A_Imp V_Imp F_Imp*r_p2];
N=length(A_Imp);
B_Imp=[ones(N,1) Data_Imp(:,1:2)]\Data_Imp(:,3);
a_Imp=B_Imp(2);
b_Imp=B_Imp(3);
c_Imp=B_Imp(1);
x_Imp=[];
y_Imp=[];
z_Imp=[];
for i=-300:1:300
    for j=-10:1/30:10
        x_Imp(end+1)=i;
        y_Imp(end+1)=j;
        z_Imp(end+1)=a_Imp*i+b_Imp*j+c_Imp;
    end
end






A=table2array(Ameas)*pi/180;
F=table2array(Fmeas);
V=table2array(Vmeas)*pi/180;

Data=[A V F*r_p2];
N=length(A);
B=[ones(N,1) Data(:,1:2)]\Data(:,3);
a=B(2);
b=B(3);
c=B(1);
x=[];
y=[];
z=[];
for i=-300:1:300
    for j=-10:1/30:10
        x(end+1)=i;
        y(end+1)=j;
        z(end+1)=a*i+b*j+c;
    end
end




figure
plot3(A,V,F*r_p2,'o')
grid on
hold on
plot3(A_Imp,V_Imp,F_Imp*r_p2,'*')
xlabel('Angular Acc [rad/s^2]')
ylabel('Angular Vel [rad/s]')
zlabel('Torque [Nm]')
title('Transparency Plane I=0.000001 kg*m^2')
plot3(A_adm,V_adm,F_adm*r_p2,'*')
xlabel('Angular Acc [rad/s^2]')
ylabel('Angular Vel [rad/s]')
zlabel('Torque [Nm]')
xlim([-Asaturation Asaturation])
ylim([-Vsaturation Vsaturation])
zlim([-Tsaturation Tsaturation])
plot3(x,y,z,'-','LineWidth',1)
plot3(x_adm,y_adm,z_adm,'-','LineWidth',1)
plot3(x_Imp,y_Imp,z_Imp,'-','LineWidth',1)
legend('Measurement no control','Measurement feed-forward compensation','Measurement admittance control','Transparency plane no control',...
    'Transparency Plane Admittance control','Transparency plane feed-forward compensation',...
    'Location','NorthOutside')

maxFadm=max(F_adm*r_p2);
minFadm=min(F_adm*r_p2);
maxF=max(F*r_p2);
minF=min(F*r_p2);
time=0:1/2000:10.0005;
time=time';
time1=0:1/1000:10.001;
time1=time1';

figure

subplot(3,1,1)
plot(time,F_adm*r_p2,'k')
title('Admittance Control: Torque Velocity and Acceleration')

ylabel('\tau_{adm} [Nm]')
xlabel('time [s]')
xlim([0,10])

subplot(3,1,2)
plot(time,V_adm,'b')
ylabel('V_{adm} [rad/s]')
xlabel('time [s]')
xlim([0,10])

subplot(3,1,3)
plot(time,A_adm,'g')
ylabel('A_{adm} [rad/s^2]')
xlabel('time [s]')
xlim([0,10])

figure
subplot(3,1,1)
plot(time,F*r_p2,'k')
title('No Control: Torque Velocity and Acceleration')
ylabel('\tau [Nm]')
xlabel('time [s]')
xlim([0,10])

subplot(3,1,2)
plot(time,V,'b')
ylabel('V [rad/s]')
xlabel('time [s]')
xlim([0,10])

subplot(3,1,3)
plot(time,A,'g')
ylabel('A [rad/s^2]')
xlabel('time [s]')
xlim([0,10])

figure

subplot(3,1,1)
plot(time1,F_Imp*r_p2,'k')
title('Feed-forward compensation: Torque Velocity and Acceleration')

ylabel('\tau_{Imp} [Nm]')
xlabel('time [s]')
xlim([0,10])

subplot(3,1,2)
plot(time1,V_Imp,'b')
ylabel('V_{Imp} [rad/s]')
xlabel('time [s]')
xlim([0,10])

subplot(3,1,3)
plot(time1,A_Imp,'g')
ylabel('A_{Imp} [rad/s^2]')
xlabel('time [s]')
xlim([0,10])

%% Transparency planes with fitting toolbox
Asat=readtable('Accsat.xlsx');
Fsat=readtable('Fsat.xlsx');
Vsat=readtable('Vsat.xlsx');


A_sat=table2array(Asat)*pi/180;
F_sat=table2array(Fsat);
V_sat=table2array(Vsat)*pi/180;

Asaturation=max(abs(A_sat));
Fsaturation=max(abs(F_sat));
Vsaturation=max(abs(V_sat));
Tsaturation=Fsaturation*r_p2;


Ameas=readtable('AccMeas.xlsx');
Fmeas=readtable('Fmeas.xlsx');
Vmeas=readtable('Vmeas.xlsx');

Ameas_adm=readtable('AccMeasadm.xlsx');
Fmeas_adm=readtable('Fmeasadm.xlsx');
Vmeas_adm=readtable('Vmeasadm.xlsx');


A_adm=table2array(Ameas_adm)*pi/180;
F_adm=table2array(Fmeas_adm);
V_adm=table2array(Vmeas_adm)*pi/180;

A=table2array(Ameas)*pi/180;
F=table2array(Fmeas);
V=table2array(Vmeas)*pi/180;

[FitImp,gofimp]=fit([A_Imp V_Imp],F_Imp*r_p2,'poly11');
[FitAdm,gofadm]=fit([A_adm V_adm],F_adm*r_p2,'poly11');
[FitNoControl,gof]=fit([A V],F*r_p2,'poly11');
figure
plot(FitAdm)
hold on
plot(FitNoControl)
plot(FitImp)
plot3(A,V,F*r_p2,'o')
plot3(A_adm,V_adm,F_adm*r_p2,'*')
plot3(A_Imp,V_Imp,F_Imp*r_p2,'*')
xlabel('Angular Acc [rad/s^2]')
ylabel('Angular Vel [rad/s]')
zlabel('Torque [Nm]')
title('Transparency Plane I=0.000001 kg*m^2')
xlim([-Asaturation Asaturation])
ylim([-Vsaturation Vsaturation])
zlim([-Tsaturation Tsaturation])




%% Saturation Force Vel Acc for low admittance case

Asat=readtable('Accsat.xlsx');
Fsat=readtable('Fsat.xlsx');
Vsat=readtable('Vsat.xlsx');


A_sat=table2array(Asat)*pi/180;
F_sat=table2array(Fsat);
V_sat=table2array(Vsat)*pi/180;

Asaturation=max(abs(A_sat));
Fsaturation=max(abs(F_sat));
Vsaturation=max(abs(V_sat));
Tsaturation=Fsaturation*r_p2;

%% Low admittance (almost meaningless)

Ameas=readtable('AccMeas.xlsx');
Fmeas=readtable('Fmeas.xlsx');
Vmeas=readtable('Vmeas.xlsx');

Ameas_adm=readtable('AccMeasadm.xlsx');
Fmeas_adm=readtable('Fmeasadm.xlsx');
Vmeas_adm=readtable('Vmeasadm.xlsx');


A_adm=table2array(Ameas_adm)*pi/180;
F_adm=table2array(Fmeas_adm);
V_adm=table2array(Vmeas_adm)*pi/180;

Data_adm=[A_adm V_adm F_adm*r_p2];
N=length(A_adm);
B_adm=[ones(N,1) Data_adm(:,1:2)]\Data_adm(:,3);
a_adm=B_adm(2);
b_adm=B_adm(3);
c_adm=B_adm(1);
x_adm=[];
y_adm=[];
z_adm=[];
for i=-300:1:300
    for j=-10:1/30:10
        x_adm(end+1)=i;
        y_adm(end+1)=j;
        z_adm(end+1)=a_adm*i+b_adm*j+c_adm;
    end
end

A=table2array(Ameas)*pi/180;
F=table2array(Fmeas);
V=table2array(Vmeas)*pi/180;

Data=[A V F*r_p2];
N=length(A);
B=[ones(N,1) Data(:,1:2)]\Data(:,3);
a=B(2);
b=B(3);
c=B(1);
x=[];
y=[];
z=[];
for i=-300:1:300
    for j=-10:1/30:10
        x(end+1)=i;
        y(end+1)=j;
        z(end+1)=a*i+b*j+c;
    end
end




figure
plot3(A,V,F*r_p2,'o')
hold on
xlabel('Angular Acc [rad/s^2]')
ylabel('Angular Vel [rad/s]')
zlabel('Torque [Nm]')
title('Transparency Plane I=0.5 kg*m^2')
plot3(A_adm,V_adm,F_adm*r_p2,'*')
xlabel('Angular Acc [rad/s^2]')
ylabel('Angular Vel [rad/s]')
zlabel('Torque [Nm]')
xlim([-Asaturation Asaturation])
ylim([-Vsaturation Vsaturation])
zlim([-Tsaturation Tsaturation])
plot3(x,y,z,'-','LineWidth',1)
plot3(x_adm,y_adm,z_adm,'-','LineWidth',1)
legend('Measurement no control','Measurement admittance control','Transparency plane no control',...
    'Transparency Plane Admittance control','Location','NorthOutside')
%legend('Measurement no control','Transparency plane no control','Location','NorthOutside')
%legend('Measurement admittance control','Transparency plane admittance control','Location','NorthOutside')

%% MBK-Plot

figure
title('MBK-Plot')

xlabel('Virtual Damping B [Nms/rad]')
ylabel('Virtual Stiffness K [Nm/rad]')
grid on
hold on
%Measurement 1
I1 = 0.001;
B1 = [0     0.001   0.005   0.01    0.012   0.015   0.02    0.025   0.027   0.0285];
K1 = [0.00005     0.19    0.39    0.42    0.38    0.32    0.22    0.1     0.05    0];

I2 = 0.01;
B2 = [0     0.001   0.005   0.01    0.05    0.1     0.2     0.25    0.27   ];
K2 = [0.06    0.15    0.95    1.6     4.1     4.1     2.2     0.8     0];


I3 = 0.1;
B3 = [0    0.001    0.05    0.1     0.3     0.5     0.8   0.9   1.1     1.25    1.5     2.0     2.5     2.6 2.7 2.72];
K3 = [0.06    0.1      8.5     16      30      35      41    45    41        39      34     26      10     5   1   0];

plot(B1, K1, 'b-o', 'LineWidth',2.0)
plot(B2, K2, 'k-o', 'LineWidth',2.0 )
plot(B3, K3, 'r-o', 'LineWidth',2.0 )
%legend('I = 0.001 kgm^2')
%legend('I = 0.001 kgm^2', 'I = 0.01 kgm^2')
legend('I = 0.001 kgm^2', 'I = 0.01 kgm^2', 'I = 0.1 kgm^2')
%% Zwidth calculation
N=length(B1);
Z_width1=0;
for i=1:N-1
    Z_width1=Z_width1+((K1(i)+K1(i+1))*(B1(i+1)-B1(i)))/2; % [N^2*m^2*s/rad^2]
end

N=length(B2);
Z_width2=0;
for i=1:N-1
    Z_width2=Z_width2+((K2(i)+K2(i+1))*(B2(i+1)-B2(i)))/2; 
end

N=length(B3);
Z_width3=0;
for i=1:N-1
    Z_width3=Z_width3+((K3(i)+K3(i+1))*(B3(i+1)-B3(i)))/2; 
end




















