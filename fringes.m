clc,clear,close all
tic
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;
%all /m
d = 1.92;%satellite distance
L = 550;%distance
%H = 40;%figure size LARGE
H = 1;%figure size SMALL
% THE DISTANCE BETWEEN FRINGES = lambda*sqrt(L^2-d^2/4)/d
delta = lambda*sqrt(L^2-d^2/4)/d;
angle = -360:1:360;
%satellite coordinates
x01 = -d/2+delta/2;
y01 = 0;
z01 = sqrt(L^2-x01^2);

x02 = d/2+delta/2;
y02 = 0;
z02 = sqrt(L^2-x02^2);
%aim = (x01+x02)/2;
aim = 0;

step = H/300;
y = -H/2:step:H/2;
x = -H/2:step:H/2;
theta_aim1 = atand(0/(0-x01));
phi_aim1 = atand(sqrt((aim-x01)^2+0^2)/z01);
theta_aim2 = atand(0/(0-x02));
phi_aim2 = atand(sqrt((aim-x02)^2+0^2)/z02);

gain1 = beamforming2d(4,4,theta_aim1,phi_aim1,frequency,angle);
gain2 = beamforming2d(4,4,theta_aim2,phi_aim2,frequency,angle);
gain1 = gain1.';
gain2 = gain2.';
gain1 = abs(gain1)./max(abs(gain1));
gain2 = abs(gain2)./max(abs(gain2));
%gain每一列对应俯仰角变化 每一行对应方位角变化


%%
[X Y] = meshgrid(x,y);%the z=0 plane
cor = zeros(length(x),length(y)); 
phi_re = zeros(length(x),length(y));
wavepath_dif = zeros(length(x),length(y)); %in wavenumber
for m = 1:length(x)
    for n = 1:length(y)
        R1 = sqrt((y(n)-y01).^2+(x(m)-x01).^2+z01.^2);
        R2 = sqrt((y(n)-y02).^2+(x(m)-x02).^2+z02.^2);
        wavepath_dif(m,n) =  (R1-R2)/lambda;
        theta_all1 = atand(y(n)/(x(m)-x01));
        phi_all1 = atand(((x(m)-x01)^2+y(n)^2)/z01);
        phi_re(m,n) = phi_all1;
        theta_all2 = atand(y(n)/(x(m)-x02));
        phi_all2 = atand(sqrt((x(m)-x02)^2+y(n)^2)/z01);
        [aa,idx_phi1] = min(abs(phi_all1-angle));
        [bb,idx_theta1] = min(abs(theta_all1-angle));
        [cc,idx_phi2] = min(abs(phi_all2-angle));
        [dd,idx_theta2] = min(abs(theta_all2-angle));
        A1 = gain1(idx_phi1,idx_theta1).*exp(1i*2*pi/lambda*R1);
        A2 = gain2(idx_phi2,idx_theta2).*exp(1i*2*pi/lambda*R2);
        I = (A1+A2)*conj(A1+A2);
        cor(m,n) = I;
    end
end
cor = 10*log10(cor);
%cor(cor<0)=0;
%%
figure
meshc(X,Y,cor')
title('power distribution/db'); 
xlabel('x-axis','fontsize',14);  
ylabel('y-axis','fontsize',14); 
zlabel('power gain/db','fontsize',14)

figure 
plot(x,cor(:,1))
title('slice of power distribution/db'); %x_axis
xlabel('x-axis','fontsize',14);  
ylabel('power gain/db','fontsize',14)
toc
%%
figure
meshc(X,Y,phi_re)
title('phi distribution/db'); 
xlabel('x-axis','fontsize',14);  
ylabel('y-axis','fontsize',14); 
zlabel('\phi/degree','fontsize',14)

figure
plot(x,wavepath_dif(:,1))
title('wavenumber distribution/db'); 
xlabel('x-axis','fontsize',14);   
%%
a= find(mean(cor,2)>2);
b=x(a);
plot(1:length(b),b,'o')
