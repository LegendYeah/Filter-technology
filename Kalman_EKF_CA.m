%% kalman filter 
clear
Sat1_pos_ang=importdata('Sat1_pos_ang.csv') % Time (EpSec)	x (km)	y (km)	z (km)	Angle (deg)	Angle (deg)
Sat2_pos_ang=importdata('Sat2_pos_ang.csv') % Time (EpSec)	x (km)	y (km)	z (km)	Angle (deg)	Angle (deg)
Tar_pos_vel=importdata('Tar_pos_vel.csv')   % Time (EpSec)	x (km)	y (km)	z (km)	vx (km/sec)	vy (km/sec)	vz (km/sec)

Sat1_pos_ang=Sat1_pos_ang.data;% Measure value1
Sat2_pos_ang=Sat2_pos_ang.data;% Measure value1
Tar_pos_vel=Tar_pos_vel.data; % True value
num=size(Sat1_pos_ang,1);% number of the data 
delta_T=0.1;
%% filter demo 2 (CA)
% Initialization
Xekf_a0=1*[1;1;1];
% Xekf(:,1)=; % initial value for filter,column vector
Tar=[4180.411;-3898.840;2819.94;-2.3034;-3.3771;6.2815];
% Xekf(:,1)=[Tar_pos_vel(1,2:7)';Xekf_a0];
Xekf(:,1)=[Tar;Xekf_a0];
% Xekf(:,1)=[P1(1,:)';V';Xekf_a0];
Phi=[eye(3),delta_T*(eye(3)),delta_T^2/2*(eye(3));zeros(3),eye(3),delta_T*(eye(3));zeros(3),zeros(3),eye(3)];
Gamma=[eye(3);zeros(3);zeros(3)];
% Q=1e-3*eye(9);
Q=[5e-3*eye(3) zeros(3) zeros(3);zeros(3) zeros(3) zeros(3);zeros(3) zeros(3) 1e-3*eye(3)]% km system noise matrix[5 0 0; 0 5 0; 0 0 5]
% R=1/57.3*[0.02 0 0 0; 0 0.02 0 0; 0 0 0.02 0; 0 0 0 0.02]; % degree observe noise
R=1e-2*[0.02 0 0 0; 0 0.02 0 0; 0 0 0.02 0; 0 0 0 0.02];
P0=1*eye(9);% covariance matrix
time(1)=0;
for i=2:num % kalman filter process
Xn=Phi*Xekf(:,i-1); % one step prediction
P1=Phi*P0*Phi'+Q; %covariance matrix prediction
% utilizing one step prediction to solve H
x_sat1=Sat1_pos_ang(i,2);y_sat1=Sat1_pos_ang(i,3); z_sat1=Sat1_pos_ang(i,4);
x_sat2=Sat2_pos_ang(i,2);y_sat2=Sat2_pos_ang(i,3); z_sat2=Sat2_pos_ang(i,4);
x1=Xn(1)-x_sat1; y1=Xn(2)-y_sat1; z1=Xn(3)-z_sat1;
x2=Xn(1)-x_sat2; y2=Xn(2)-y_sat2; z2=Xn(3)-z_sat2;
% Jacobi matrix
H=[-y1/(x1^2+y1^2) x1/(x1^2+y1^2)  0 0 0 0 0 0 0;
    -x1*z1/(sqrt(x1^2+y1^2)*(x1^2+y1^2+z1^2)) -y1*z1/(sqrt(x1^2+y1^2)*(x1^2+y1^2+z1^2)) sqrt(x1^2+y1^2)/(x1^2+y1^2+z1^2) 0 0 0 0 0 0;
   -y2/(x2^2+y2^2) x2/(x2^2+y2^2)  0 0 0 0 0 0 0;
    -x2*z2/(sqrt(x2^2+y2^2)*(x2^2+y2^2+z2^2)) -y2*z2/(sqrt(x2^2+y2^2)*(x2^2+y2^2+z2^2)) sqrt(x2^2+y2^2)/(x2^2+y2^2+z2^2) 0 0 0 0 0 0];
K=P1*H'*pinv(H*P1*H'+R);% gain matrix
Z=[Sat1_pos_ang(i,5:6),Sat2_pos_ang(i,5:6)]'/(180/pi);%True observation£¬ convert ¡¾deg¡¿into ¡¾rad¡¿
Z_pre=[atan2(y1,x1),atan(z1/sqrt(x1^2+y1^2)),atan2(y2,x2),atan(z2/sqrt(x2^2+y2^2))]';% therotical observe value
Xekf(:,i)=Xn+K*(Z-Z_pre); % filter value
P0=(eye(9)-K*H)*P1;% simplify formula
observe_error(:,i-1)=Z-Z_pre;
time(i)=(i-1)*delta_T;
end
% filter error
error_position=Tar_pos_vel(:,2:4)'-Xekf(1:3,:);
error_velocity=Tar_pos_vel(:,5:7)'-Xekf(4:6,:);

% figrue plot
figure(1)
plot3(Tar_pos_vel(:,2),Tar_pos_vel(:,3),Tar_pos_vel(:,4),'r','LineWidth',2) %target trajectory
hold on 
plot3(Xekf(1,:),Xekf(2,:),Xekf(3,:),'b','LineWidth',2) 
legend('Target true trajectory','Target filter trajectory')
grid on
xlabel('x/km','FontSize',12,'FontWeight','normal')
ylabel('y/km','FontSize',12,'FontWeight','normal')
zlabel('z/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
figure(2) % error
subplot(3,1,1)
plot(time,error_position(1,:),'k','LineWidth',1.5)
grid on 
title('x')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,error_position(2,:),'k','LineWidth',1.5)
grid on 
title('y')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,error_position(3,:),'k','LineWidth',1.5)
grid on 
title('z')
xlabel('t/s','FontSize',12,'FontWeight','normal')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
figure(3) % error
subplot(3,1,1)
plot(time,error_velocity(1,:),'k','LineWidth',1.5)
grid on 
title('vx')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,error_velocity(2,:),'k','LineWidth',1.5)
grid on 
title('vy')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,error_velocity(3,:),'k','LineWidth',1.5)
grid on 
title('vz')
xlabel('t/s','FontSize',12,'FontWeight','normal')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');

figure(4) % acceleration
subplot(3,1,1)
plot(time,Xekf(7,:),'k','LineWidth',1.5)
grid on 
title('ax')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,Xekf(8,:),'k','LineWidth',1.5)
grid on 
title('ay')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,Xekf(9,:),'k','LineWidth',1.5)
grid on 
title('az')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
xlabel('t/s','FontSize',12,'FontWeight','normal')


