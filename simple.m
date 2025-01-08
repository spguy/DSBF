clc,clear,close all
tic
c = 3e8;
frequency = 13.7e9;
lambda = c/frequency;%波长
%all /m
d =3/137;%satellite distance 100m
D = 0.5;%distance 100km
delta = lambda*D/d;
angle = -360:1:360;
%satellite coordinates
x01 = d;
y01 = d;
z01 = D;

x02 = -d;
y02 = d;
z02 = D;

x03 = d;
y03 = -d;
z03 = D;

x04 = -d;
y04 = -d;
z04 = D;

aim = 0;

H = 200;
step = H/2001;
y = -H/2:step:H/2;
x = -H/2:step:H/2;
phi1 = -0.36*pi;
phi2 = 0.98*pi;
phi3 = -0.71*pi;
phi4 = 0.11*pi;
A1 = 1;
A2 = 1;
A3 = 0;
A4 = 0;

% phi1 = 0;
% phi2 = 0;
% phi3 = 0;
% phi4 = 0;
% A1 = 1;
% A2 = 1;
% A3 = 1;
% A4 = 1;

t = 0:1e-11:1e-9;
%产生两个正交的线极化波，每个线极化波由x方向和y方向的同相位阵子合成而得
A_x = 1;
A_y = 1;
s1_x = A_x*exp(1i*(2*pi*frequency*t+phi1));
s1_y = A_y*exp(1i*(2*pi*frequency*t+phi1));
s2_x = -A_y*exp(1i*(2*pi*frequency*t+phi2));
s2_y = A_x*exp(1i*(2*pi*frequency*t+phi2));

s1 = A1*exp(1i*(2*pi*frequency*t+phi1));
s2 = A2*exp(1i*(2*pi*frequency*t+phi2));
s3 = A3*exp(1i*(2*pi*frequency*t+phi3));
s4 = A4*exp(1i*(2*pi*frequency*t+phi4));
cor = zeros(length(x),length(y));
cor_x = zeros(length(x),length(y));
cor_y = zeros(length(x),length(y));
signal_map = zeros(length(t),length(x),length(y));
for m = 1:length(x)
    for n = 1:length(y)
        R1 = sqrt((y(n)-y01).^2+(x(m)-x01).^2+z01.^2);
        R2 = sqrt((y(n)-y02).^2+(x(m)-x02).^2+z02.^2);
        R3 = sqrt((y(n)-y03).^2+(x(m)-x03).^2+z03.^2);
        R4 = sqrt((y(n)-y04).^2+(x(m)-x04).^2+z04.^2);
        delta_phi1 = R1*2*pi/lambda;     
        delta_phi2 = R2*2*pi/lambda; 
        delta_phi3 = R3*2*pi/lambda; 
        delta_phi4 = R4*2*pi/lambda; 
        signal_receive = (s1)*exp(1i*delta_phi1)+(s2)*exp(1i*delta_phi2)+(s3)*exp(1i*delta_phi3)+(s4)*exp(1i*delta_phi4);
        %signal_receive_x = s1_x*exp(1i*delta_phi1)+s2_x*exp(1i*delta_phi2);
        %signal_receive_y = s1_y*exp(1i*delta_phi1)+s2_y*exp(1i*delta_phi2);
        cor(m,n) = (signal_receive*signal_receive'/length(t));
        %cor_x(m,n) = (signal_receive_x*signal_receive_x'/length(t));
        %cor_y(m,n) = (signal_receive_y*signal_receive_y'/length(t));
    end
end
%%
cor(cor<0) = 0;
cor = 10*log10(cor);
cor(cor<-10)=-10;
%%
figure
[X Y] = meshgrid(x,y); 
meshc(X,Y,cor')
xlabel('x-axis(m)','fontsize',64,'FontName','Times New Roman');  
ylabel('y-axis(m)','fontsize',64,'FontName','Times New Roman'); 
zlabel('Power Gain/dB','fontsize',64,'FontName','Times New Roman')
zlim([-10,9])
%%
node_x = randperm(length(x),10);
node_y = randperm(length(y),10); 
figure
imagesc(X(1, :), Y(:, 1), cor);
xlabel('x-axis(m)','fontsize',48,'FontName','Times New Roman');  
ylabel('y-axis(m)','fontsize',48,'FontName','Times New Roman'); 
%%
hold on
plot(x(node_x), y(node_y), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
