%% kalman filter  Singer
% load data 
clear
Sat1_pos_ang=importdata('Sat1_pos_ang.csv') % Time (EpSec)	x (km)	y (km)	z (km)	Angle (deg)	Angle (deg)
Sat2_pos_ang=importdata('Sat2_pos_ang.csv') % Time (EpSec)	x (km)	y (km)	z (km)	Angle (deg)	Angle (deg)
Tar_pos_vel=importdata('Tar_pos_vel.csv')   % Time (EpSec)	x (km)	y (km)	z (km)	vx (km/sec)	vy (km/sec)	vz (km/sec)

Sat1_pos_ang=Sat1_pos_ang.data;% Measure value1
Sat2_pos_ang=Sat2_pos_ang.data;% Measure value1
Tar_pos_vel=Tar_pos_vel.data; % True value
num=size(Sat1_pos_ang,1);% number of the data 
delta_T=0.1;
% choose different motion frequency as consideration
alpha_set=[1;1/20;1/60]
for nn=1:3
alpha=alpha_set(nn);
A0=exp(-1*alpha*delta_T);
A=(1-exp(-1*alpha*delta_T))/alpha;
A2=(alpha*delta_T-1+exp(-1*alpha*delta_T))/alpha^2;
% covariance parameters
q11=(1+2*alpha*delta_T-2*alpha^2*delta_T^2+2/3*alpha^3*delta_T^3-...
    exp(-2*alpha*delta_T)-4*alpha*delta_T*exp(-alpha*delta_T))/(2*alpha^5);
q12=(1-2*alpha*delta_T+alpha^2*delta_T^2-2*exp(-alpha*delta_T)+exp(-2*alpha*delta_T)+2*alpha*delta_T*exp(-alpha*delta_T))/(2*alpha^4);
q21=q12;
q13=(1-exp(-2*alpha*delta_T)-2*alpha*delta_T*exp(-alpha*delta_T))/(2*alpha^3);
q31=q13;
q22=(-3+2*alpha*delta_T-exp(-2*alpha*delta_T)+4*exp(-alpha*delta_T))/(2*alpha^3);
q23=(1-2*exp(-alpha*delta_T)+exp(-2*alpha*delta_T))/(2*alpha^3);
q32=q23;
q33=(1-exp(-2*alpha*delta_T))/(2*alpha);
%% filter demo 3 (Singer model)
% Initialization
Xekf(:,1)=[4180.411;-2.3034;1;-3898.840;-3.3771;1;2819.94;6.2815;1]; % intial value x vx ax y vy ay z vz az
F=[1 delta_T A2;0 1 A;0 0 A0];
Phi=[F zeros(3) zeros(3);zeros(3) F,zeros(3);zeros(3),zeros(3) F];
sigma=1;
qq=2*alpha*sigma^2*[q11 q12 q13;q21 q22 q23;q31 q32 q33] % km system noise matrix
Q=[qq zeros(3) zeros(3);zeros(3) qq,zeros(3);zeros(3),zeros(3) qq];
R=1e-2*[0.02 0 0 0; 0 0.02 0 0; 0 0 0.02 0; 0 0 0 0.02]; % degree observe noise
P0=eye(9);% covariance matrix
for i=2:num % kalman filter process
Xn=Phi*Xekf(:,i-1); % one step prediction
P1=Phi*P0*Phi'+Q; %covariance matrix prediction
% utilizing one step prediction to solve H
x_sat1=Sat1_pos_ang(i,2);y_sat1=Sat1_pos_ang(i,3); z_sat1=Sat1_pos_ang(i,4);
x_sat2=Sat2_pos_ang(i,2);y_sat2=Sat2_pos_ang(i,3); z_sat2=Sat2_pos_ang(i,4);
x1=Xn(1)-x_sat1; y1=Xn(4)-y_sat1; z1=Xn(7)-z_sat1;
x2=Xn(1)-x_sat2; y2=Xn(4)-y_sat2; z2=Xn(7)-z_sat2;
% Jacobi matrix
H=[-y1/(x1^2+y1^2) 0 0 x1/(x1^2+y1^2) 0 0 0 0 0  ;
    -x1*z1/(sqrt(x1^2+y1^2)*(x1^2+y1^2+z1^2)) 0 0 -y1*z1/(sqrt(x1^2+y1^2)*(x1^2+y1^2+z1^2)) 0 0 sqrt(x1^2+y1^2)/(x1^2+y1^2+z1^2) 0 0 ;
   -y2/(x2^2+y2^2) 0 0 x2/(x2^2+y2^2) 0 0  0 0 0 ;
    -x2*z2/(sqrt(x2^2+y2^2)*(x2^2+y2^2+z2^2)) 0 0 -y2*z2/(sqrt(x2^2+y2^2)*(x2^2+y2^2+z2^2)) 0 0 sqrt(x2^2+y2^2)/(x2^2+y2^2+z2^2) 0 0 ];
K=P1*H'*inv(H*P1*H'+R);% gain matrix
Z=[Sat1_pos_ang(i,5:6),Sat2_pos_ang(i,5:6)]'/(180/pi);%True observation£¬ convert ¡¾deg¡¿into ¡¾rad¡¿
Z_pre=[atan2(y1,x1),atan(z1/sqrt(x1^2+y1^2)),atan2(y2,x2),atan(z2/sqrt(x2^2+y2^2))]';% therotical observe value
Xekf(:,i)=Xn+K*(Z-Z_pre); % filter value
P0=(eye(9)-K*H)*P1;% covariance matrix
observe_error(:,i-1)=Z-Z_pre;
time(i)=(i-1)*delta_T;
end
% filter error& accleration
x(nn,:)=Xekf(1,:);
y(nn,:)=Xekf(4,:);
z(nn,:)=Xekf(7,:);
error_position=[Tar_pos_vel(:,2)'-Xekf(1,:);Tar_pos_vel(:,3)'-Xekf(4,:);Tar_pos_vel(:,4)'-Xekf(7,:)];
error_x(nn,:)=error_position(1,:); % i=1 alpha=1,i=2 alpha=1/20; i=3 alpha=1/60;
error_y(nn,:)=error_position(2,:);
error_z(nn,:)=error_position(3,:);
error_velocity=[Tar_pos_vel(:,5)'-Xekf(2,:);Tar_pos_vel(:,6)'-Xekf(5,:);Tar_pos_vel(:,7)'-Xekf(8,:)];
error_vx(nn,:)=error_velocity(1,:); % i=1 alpha=1,i=2 alpha=1/20; i=3 alpha=1/60;
error_vy(nn,:)=error_velocity(2,:);
error_vz(nn,:)=error_velocity(3,:);
ax(nn,:)=Xekf(3,:);
ay(nn,:)=Xekf(6,:);
az(nn,:)=Xekf(9,:);
end
error_pos_singer(1)=norm([mean(error_x(1,2000:7000)),mean(error_y(1,2000:7000)),mean(error_z(1,2000:7000))]);
error_pos_singer(2)=norm([mean(error_x(2,2000:7000)),mean(error_y(2,2000:7000)),mean(error_z(2,2000:7000))]);
error_pos_singer(3)=norm([mean(error_x(3,2000:7000)),mean(error_y(3,2000:7000)),mean(error_z(3,2000:7000))]);
% figrue plot
figure(1) % trajectory plot
plot3(Tar_pos_vel(:,2),Tar_pos_vel(:,3),Tar_pos_vel(:,4),'r','LineWidth',2) %target trajectory
hold on 
plot3(x(1,:),y(1,:),z(1,:),'b:','LineWidth',2) 
hold on 
plot3(x(2,:),y(2,:),z(2,:),'b--','LineWidth',2) 
hold on
plot3(x(3,:),y(3,:),z(3,:),'b','LineWidth',2) 
legend('Target true trajectory','Target filter trajectory \alpha=1','Target filter trajectory \alpha=1/20','Target filter trajectory \alpha=1/60')
grid on
xlabel('x/km','FontSize',12,'FontWeight','normal')
ylabel('y/km','FontSize',12,'FontWeight','normal')
zlabel('z/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
figure(2) % error position
subplot(3,1,1)
plot(time,error_x(1,:),'k','LineWidth',1.2)
hold on
plot(time,error_x(2,:),'k:','LineWidth',1.5)
hold on
plot(time,error_x(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('x')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,error_y(1,:),'k','LineWidth',1.2)
hold on
plot(time,error_y(2,:),'k:','LineWidth',1.5)
hold on
plot(time,error_y(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('y')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,error_z(1,:),'k','LineWidth',1.2)
hold on
plot(time,error_z(2,:),'k:','LineWidth',1.5)
hold on
plot(time,error_z(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('z')
xlabel('t/s','FontSize',12,'FontWeight','normal')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
figure(3) % error velocity 
subplot(3,1,1)
plot(time,error_vx(1,:),'k','LineWidth',1.2)
hold on
plot(time,error_vx(2,:),'k:','LineWidth',1.5)
hold on
plot(time,error_vx(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('vx')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,error_vy(1,:),'k','LineWidth',1.2)
hold on
plot(time,error_vy(2,:),'k:','LineWidth',1.5)
hold on
plot(time,error_vy(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('vy')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,error_vz(1,:),'k','LineWidth',1.2)
hold on
plot(time,error_vz(2,:),'k:','LineWidth',1.5)
hold on
plot(time,error_vz(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('vz')
xlabel('t/s','FontSize',12,'FontWeight','normal')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');

figure(4) % target acceleration 
subplot(3,1,1)
plot(time,ax(1,:),'k','LineWidth',1.2)
hold on
plot(time,ax(2,:),'k:','LineWidth',1.5)
hold on
plot(time,ax(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('ax')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,ay(1,:),'k','LineWidth',1.2)
hold on
plot(time,ay(2,:),'k:','LineWidth',1.5)
hold on
plot(time,ay(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('ay')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,az(1,:),'k','LineWidth',1.2)
hold on
plot(time,az(2,:),'k:','LineWidth',1.5)
hold on
plot(time,az(3,:),'k--','LineWidth',1.5)
grid on 
legend('\alpha=1','\alpha=1/20','\alpha=1/60')
title('az')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
xlabel('t/s','FontSize',12,'FontWeight','normal')
