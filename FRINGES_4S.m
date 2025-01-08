clc,clear,close all
tic
%satellite orbit is the x axis
c = 3e8;
frequency = 3.5e8;
lambda = c/frequency;%wavelength
%all /km
d = 1.92;%satellite distance
L = 55;%orbit radius
beamwidth = 10;%degree
% calculation of covering area on the ground
%H = L*tand(beamwidth/2);%figure size
H = 100;
step = H/2001;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;
%lambda = lambda/1000; %change into /km
angle = -360:1:360;
% THE DISTANCE BETWEEN FRINGES = lambda*sqrt(L^2-d^2/4)/d
delta_compen = lambda*sqrt(L^2-d^2/4)/d;
%initial satellites locations
x01 = d/2;
y01 = 0;
z01 = sqrt(L^2-x01^2-y01^2);
  
x02 = -d/2;
y02 = 0;
z02 = sqrt(L^2-x02^2-y02^2);
%aim1 = 0;
aim1 = 0;

x03 = 0;
y03 = d/2;
z03 = sqrt(L^2-x03^2-y03^2);

x04 = 0;
y04 = -d/2;
z04 = sqrt(L^2-x04^2-y04^2);
%aim2 = 0;
aim2 = 0;

step = H/2001;
y = -H:step:H;
x = -H:step:H;

%initial beam angles for satellites arrays 
theta_aim1 = atand(0/(aim1-x01));
phi_aim1 = atand(sqrt((aim1-x01)^2+0^2)/z01);

theta_aim2 = atand(0/(aim1-x02));
phi_aim2 = atand(sqrt((aim1-x02)^2+0^2)/z02);

theta_aim3 = atand(0/(aim2-x03));
phi_aim3 = atand(sqrt((aim2-x03)^2+0^2)/z03);

theta_aim4 = atand(0/(0-x04));
phi_aim4 = atand(sqrt((aim2-x04)^2+0^2)/z04);

gain1 = beamforming2d(8,8,theta_aim1,phi_aim1,frequency,angle);
gain2 = beamforming2d(8,8,theta_aim2,phi_aim2,frequency,angle);
gain3 = beamforming2d(8,8,theta_aim3,phi_aim3,frequency,angle);
gain4 = beamforming2d(8,8,theta_aim4,phi_aim4,frequency,angle);
gain1 = gain1.';
gain2 = gain2.';
gain3 = gain3.';
gain4 = gain4.';
gain1 = abs(gain1)./max(abs(gain1));
gain2 = abs(gain2)./max(abs(gain2));
gain3 = abs(gain3)./max(abs(gain3));
gain4 = abs(gain4)./max(abs(gain4));
%gain每一列对应俯仰角变化 每一行对应方位角变化
% 波束在地面的覆盖半径 = sin(波束宽度/2)*L

cor = zeros(length(x),length(y)); 
for m = 1:length(x)
    for n = 1:length(y)
        R1 = sqrt((y(n)-y01).^2+(x(m)-x01).^2+z01.^2);
        R2 = sqrt((y(n)-y02).^2+(x(m)-x02).^2+z02.^2);
        R3 = sqrt((y(n)-y03).^2+(x(m)-x03).^2+z03.^2);
        R4 = sqrt((y(n)-y04).^2+(x(m)-x04).^2+z04.^2);
        theta_all1 = atand(y(n)/(x(m)-x01));
        phi_all1 = atand(((x(m)-x01)^2+y(n)^2)/z01);
        theta_all2 = atand(y(n)/(x(m)-x02));
        phi_all2 = atand(sqrt((x(m)-x02)^2+y(n)^2)/z02);
        theta_all3 = atand(y(n)/(x(m)-x03));
        phi_all3 = atand(((x(m)-x03)^2+y(n)^2)/z03);
        theta_all4 = atand(y(n)/(x(m)-x04));
        phi_all4 = atand(sqrt((x(m)-x04)^2+y(n)^2)/z04);
        [aa,idx_phi1] = min(abs(phi_all1-angle));
        [bb,idx_theta1] = min(abs(theta_all1-angle));
        [cc,idx_phi2] = min(abs(phi_all2-angle));
        [dd,idx_theta2] = min(abs(theta_all2-angle));
        [aa,idx_phi3] = min(abs(phi_all3-angle));
        [bb,idx_theta3] = min(abs(theta_all3-angle));
        [cc,idx_phi4] = min(abs(phi_all4-angle));
        [dd,idx_theta4] = min(abs(theta_all4-angle));
        A1 = gain1(idx_phi1,idx_theta1).*exp(1i*2*pi/lambda*R1);
        A2 = gain2(idx_phi2,idx_theta2).*exp(1i*2*pi/lambda*R2);
        A3 = gain1(idx_phi3,idx_theta3).*exp(1i*2*pi/lambda*R3);
        A4 = gain2(idx_phi4,idx_theta4).*exp(1i*2*pi/lambda*R4);
        I = (A1+A2+A3)*conj(A1+A2+A3);
        cor(m,n) = I;
    end
end
node_x = [];
node_y = [];
[score,cor] = power_calculation(sat_num,sat_loc_x,sat_loc_y,phi,A,x,y,H,node_x,node_y,D)
powerdistribution(cor,x,y,node_x,node_y,x_ground,y_ground)

