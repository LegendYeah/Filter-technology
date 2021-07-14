%% 求空间中两点的交点坐标，作为目标位置初值
% load data 
clear
Sat1_pos_ang=importdata('Sat1_pos_ang.csv') % Time (EpSec)	x (km)	y (km)	z (km)	Angle (deg)	Angle (deg)
Sat2_pos_ang=importdata('Sat2_pos_ang.csv') % Time (EpSec)	x (km)	y (km)	z (km)	Angle (deg)	Angle (deg)
Tar_pos_vel=importdata('Tar_pos_vel.csv')   % Time (EpSec)	x (km)	y (km)	z (km)	vx (km/sec)	vy (km/sec)	vz (km/sec)

Sat1_pos_ang=Sat1_pos_ang.data;% Measure value1
Sat2_pos_ang=Sat2_pos_ang.data;% Measure value1
for i=1:2
% Sat1 initial value
initial_Sat1_Pos=Sat1_pos_ang(i,2:4);
xa=initial_Sat1_Pos(1);
ya=initial_Sat1_Pos(2);
za=initial_Sat1_Pos(3);
initial_Sat1_ang=Sat1_pos_ang(i,5:6);
eta1=initial_Sat1_ang(1);
theta1=initial_Sat1_ang(2);
% Sat2 initial value
initial_Sat2_Pos=Sat2_pos_ang(i,2:4);
xb=initial_Sat2_Pos(1);
yb=initial_Sat2_Pos(2);
zb=initial_Sat2_Pos(3);
initial_Sat2_ang=Sat2_pos_ang(i,5:6);
eta2=initial_Sat2_ang(1);
theta2=initial_Sat2_ang(2);
% Line PA
a1=-1;
b1=a1*tand(eta1);
c1=sqrt(1+(b1)^2)*tand(theta1);
% Line PB
a2=-1;
b2=a2*tand(eta2);
c2=sqrt(1+(b2)^2)*tand(theta2);
% Cross point
A=[a1 -a2;b1 -b2];
B=[xb-xa;yb-ya];
T=inv(A)*B
t1=T(1);
t2=T(2);
xp1=xa+t1*a1;
yp1=ya+t1*b1;
zp1=za+t1*c1;

xp2=xb+t2*a2;
yp2=yb+t2*b2;
zp2=zb+t2*c2;
P1(i,:)=[xp1,yp1,zp1]
P2(i,:)=[xp2,yp2,zp2]
end
% 计算V
dt=0.1;
V=(P1(2,:)-P1(1,:))/dt