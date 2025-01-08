function [score,cor] = power_calculation(sat_num,sat_loc_x,sat_loc_y,phi,x,y,H,node_x,node_y,D,frequency)
c = 3e8;
lambda = c/frequency;%波长
%satellite coordinates 
A = ones(1,sat_num);
t = 0:1e-11:1e-9;
s = zeros(sat_num,length(t));
for n = 1:sat_num
    s(n,:) = A(n)*exp(1i*(2*pi*frequency*t+phi(n)));
end
cor = zeros(length(x),length(y));
signal_map = zeros(length(t),length(x),length(y));
R = zeros(1,sat_num);
delta_phi = zeros(1,sat_num);
signal_receive = zeros(1,length(t));
%%
for m = 1:length(x)
    for n = 1:length(y)
        for k = 1:sat_num
            R(k) = sqrt((y(n)-sat_loc_y(k)).^2+(x(m)-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            signal_receive = signal_receive + s(k)*exp(1i*delta_phi(k));
        end
        cor(m,n) = signal_receive*signal_receive'/length(signal_receive);
        signal_receive = zeros(1,length(t));
    end
end
%% 计算采样点总功率
%cor(cor<0) = 0;
% signal_plot = signal_map(:,ceil(length(y) / 2),ceil(length(y) / 2));
% figure
% plot(t,signal_plot);
% title('signal at the center'); 
% legend('x','y');
% [X Y] = meshgrid(x,y); 
% figure
% meshc(X,Y,cor')
% title('power distribution/db'); 
% xlabel('x-axis','fontsize',14);  
% ylabel('y-axis','fontsize',14); 
% zlabel('power gain/db','fontsize',14)
% %%
% figure
% plot(x,cor(:,1)')
% title('power slice'); 
score = 0;
for m = 1:length(node_x)
        score = score + cor(node_x(m),node_y(m));
end
score = score/length(node_x);
end