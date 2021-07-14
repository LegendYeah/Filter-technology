%% kalman UKF VS EKF
clear
% load data 
Sat1_pos_ang=importdata('Sat1_pos_ang.csv') % Time (EpSec)	x (km)	y (km)	z (km)	Angle (deg)	Angle (deg)
Sat2_pos_ang=importdata('Sat2_pos_ang.csv') % Time (EpSec)	x (km)	y (km)	z (km)	Angle (deg)	Angle (deg)
Tar_pos_vel=importdata('Tar_pos_vel.csv')   % Time (EpSec)	x (km)	y (km)	z (km)	vx (km/sec)	vy (km/sec)	vz (km/sec)

Sat1_pos_ang=Sat1_pos_ang.data;% Measure value1
Sat2_pos_ang=Sat2_pos_ang.data;% Measure value1
Tar_pos_vel=Tar_pos_vel.data; % True value
num=size(Sat1_pos_ang,1);% number of the data 
delta_T=0.1;
%% filter  (UKF)
% parameter set
L=9; % row of status
alpha=1;
ramda=3-L;
beta=2;
%Initialization
Tar=[4180.411;-3898.840;2819.94;-2.3034;-3.3771;6.2815];
Xukf_a0=1*[1;1;1];
Xukf(:,1)=[Tar;Xukf_a0];
Phi=[eye(3),delta_T*(eye(3)),delta_T^2/2*(eye(3));zeros(3),eye(3),delta_T*(eye(3));zeros(3),zeros(3),eye(3)];
Gamma=[eye(3);zeros(3);zeros(3)];
Q=[5e-3*eye(3) zeros(3) zeros(3);zeros(3) zeros(3) zeros(3);zeros(3) zeros(3) 1e-3*eye(3)];% km system noise matrix[5 0 0; 0 5 0; 0 0 5]
% km system noise matrix
R=1e-2*[0.02 0 0 0; 0 0.02 0 0; 0 0 0.02 0; 0 0 0 0.02]; % degree observe noise
P0=eye(9);% covariance matrix
time(1)=0;
% caculate sigma weight
for j=1:2*L+1
Wm(j)=1/(2*(L+ramda));
Wc(j)=1/(2*(L+ramda));
end
Wm(1)=ramda/(L+ramda);
Wc(1)=ramda/(L+ramda)+1-alpha^2+beta;

for i=2:num % uscented kalman filter process
    xestimate=Xukf(:,i-1);
    P=P0;
    % step 1: get sigma points set
    cho=(chol(P*(L+ramda)))';
    for k=1:L
        xgamaP1(:,k)=xestimate+cho(:,k); % left sigma points set
        xgamaP2(:,k)=xestimate-cho(:,k); % right sigma points set
    end
    Xsigma=[xestimate,xgamaP1,xgamaP2]; % x sigma points
    % step 2: one step prediction
    Xsigmapre=Phi* Xsigma; % prediction: A 9*(2n+1) matrix
    % step 3: get xpre mean value and covariance
    Xpred=zeros(9,1);
    for k=1:2*L+1
        Xpred=Xpred+Wm(k)*Xsigmapre(:,k);
    end
    Ppred=zeros(9,9);
    for k=1:2*L+1
        Ppred=Ppred+Wc(k)*(Xsigmapre(:,k)-Xpred)*(Xsigmapre(:,k)-Xpred)';
    end
%   Ppred=Ppred+Gamma*Q*Gamma';
    Ppred=Ppred+Q;
    % step 4: a set of new sigma points with new UT conversion
    chor=(chol(Ppred*(L+ramda)))';
    for k=1:L
        XaugsigmaP1(:,k)=Xpred+chor(:,k);
        XaugsigmaP2(:,k)=Xpred-chor(:,k);
    end
    Xaufsigma=[Xpred XaugsigmaP1 XaugsigmaP2];% prediction: A 9*(2n+1) matrix
    % step 5: Observation prediction 
    x_sat1=Sat1_pos_ang(i,2);y_sat1=Sat1_pos_ang(i,3); z_sat1=Sat1_pos_ang(i,4);
    x_sat2=Sat2_pos_ang(i,2);y_sat2=Sat2_pos_ang(i,3); z_sat2=Sat2_pos_ang(i,4);
   for k=1:2*L+1 
       x1=Xaufsigma(1,k)-x_sat1;y1=Xaufsigma(2,k)-y_sat1;z1=Xaufsigma(3,k)-z_sat1;
       x2=Xaufsigma(1,k)-x_sat2;y2=Xaufsigma(2,k)-y_sat2;z2=Xaufsigma(3,k)-z_sat2;
       Zsigmapre(:,k)=[atan2(y1,x1),atan(z1/sqrt(x1^2+y1^2)),atan2(y2,x2),atan(z2/sqrt(x2^2+y2^2))]';
   end
  
   % step 6: 
   Zpred=zeros(4,1);
   for k=1:2*L+1 
        Zpred=Zpred+Wm(k)*Zsigmapre(:,k);
   end
   Pzz=zeros(4,4);
   for k=1:2*L+1 
   Pzz=Pzz+Wc(k)*(Zsigmapre(:,k)-Zpred)*(Zsigmapre(:,k)-Zpred)';    
   end 
   Pzz=Pzz+R; 
   Pxz=zeros(9,4);
   %step 7: Observation mean value
   for k=1:2*L+1 
   Pxz=Pxz+Wc(k)*(Xaufsigma(:,k)-Xpred)*(Zsigmapre(:,k)-Zpred)';    
   end 
   % step 8: K caculate
   K=Pxz*inv(Pzz);
   Z=[Sat1_pos_ang(i,5:6),Sat2_pos_ang(i,5:6)]'/(180/pi);
   xestimate=Xpred+K*(Z-Zpred);
   P=Ppred-K*Pzz*K';
   P0=P;
   Xukf(:,i)=xestimate;
   time(i)=(i-1)*delta_T;
end
% filter error
error_position_ukf=Tar_pos_vel(:,2:4)'-Xukf(1:3,:);
error_velocity_ukf=Tar_pos_vel(:,5:7)'-Xukf(4:6,:);
%% filter EKF
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
error_position_ekf=Tar_pos_vel(:,2:4)'-Xekf(1:3,:);
error_velocity_ekf=Tar_pos_vel(:,5:7)'-Xekf(4:6,:);
error_pos_ekf_=norm([mean(error_position_ekf(1,2000:7000)),mean(error_position_ekf(2,2000:7000)),mean(error_position_ekf(3,2000:7000))]);
error_pos_ukf_=norm([mean(error_position_ukf(1,2000:7000)),mean(error_position_ukf(2,2000:7000)),mean(error_position_ukf(3,2000:7000))]);
error_vel_ekf_=norm([mean(error_velocity_ekf(1,2000:7000)),mean(error_velocity_ekf(2,2000:7000)),mean(error_velocity_ekf(3,2000:7000))]);
error_vel_ukf_=norm([mean(error_velocity_ukf(1,2000:7000)),mean(error_velocity_ukf(2,2000:7000)),mean(error_velocity_ukf(3,2000:7000))]);
% figrue plot
figure(1)
plot3(Tar_pos_vel(:,2),Tar_pos_vel(:,3),Tar_pos_vel(:,4),'r','LineWidth',2) %target trajectory
hold on 
plot3(Xukf(1,:),Xukf(2,:),Xukf(3,:),'b','LineWidth',2) 
legend('Target true trajectory','Target filter trajectory with UKF')
grid on
xlabel('x/km','FontSize',12,'FontWeight','normal')
ylabel('y/km','FontSize',12,'FontWeight','normal')
zlabel('z/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
figure(2) % error
subplot(3,1,1)
plot(time,error_position_ekf(1,:),'bo','LineWidth',1)
hold on
plot(time,error_position_ukf(1,:),'r','LineWidth',2)
grid on 
legend('UKF','EKF')
title('x')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,error_position_ekf(2,:),'bo','LineWidth',1)
hold on
plot(time,error_position_ukf(2,:),'r','LineWidth',2)
grid on 
legend('UKF','EKF')
title('y')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,error_position_ekf(3,:),'bo','LineWidth',1)
hold on
plot(time,error_position_ukf(3,:),'r','LineWidth',2)
grid on
legend('UKF','EKF')
title('z')
xlabel('t/s','FontSize',12,'FontWeight','normal')
ylabel('distance/km','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
figure(3) % error
subplot(3,1,1)
plot(time,error_velocity_ekf(1,:),'bo','LineWidth',1)
hold on
plot(time,error_velocity_ukf(1,:),'r','LineWidth',2)
grid on 
legend('UKF','EKF')
title('vx')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,error_velocity_ekf(2,:),'bo','LineWidth',1)
hold on
plot(time,error_velocity_ukf(2,:),'r','LineWidth',2)
grid on 
legend('UKF','EKF')
title('vy')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,error_velocity_ekf(3,:),'bo','LineWidth',1)
hold on
plot(time,error_velocity_ukf(3,:),'r','LineWidth',2)
grid on 
legend('UKF','EKF')
title('vz')
xlabel('t/s','FontSize',12,'FontWeight','normal')
ylabel('velocity/km/s','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');

figure(4) % acceleration
subplot(3,1,1)
plot(time,Xekf(7,:),'bo','LineWidth',1)
hold on
plot(time,Xukf(7,:),'r','LineWidth',2)
grid on 
legend('UKF','EKF')
title('ax')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,2)
plot(time,Xekf(8,:),'bo','LineWidth',1)
hold on
plot(time,Xukf(8,:),'r','LineWidth',2)
grid on
legend('UKF','EKF')
title('ay')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
subplot(3,1,3)
plot(time,Xekf(9,:),'bo','LineWidth',1)
hold on
plot(time,Xukf(9,:),'r','LineWidth',2)
grid on 
legend('EKF','UKF')
title('az')
ylabel('\alpha/km/s^2','FontSize',12,'FontWeight','normal')
set(gca,'FontSize',12,'Fontname', 'Times New Roman','FontWeight','normal');
xlabel('t/s','FontSize',12,'FontWeight','normal')   
    