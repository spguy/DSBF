function [score] = fitness_power_nodeonly(xx,sat_num,sat_loc_x,sat_loc_y,x,y,node_x,node_y,D)
%disp(['Debug: x = ', num2str(x)]);
phi = xx(1:sat_num);
A = xx((sat_num+1):sat_num*2);
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
%satellite coordinates 
t = 0:1e-11:1e-9;
s = zeros(sat_num,length(t));
for n = 1:sat_num
    s(n,:) = A(n)*exp(1i*(2*pi*frequency*t+phi(n)));
end
cor = zeros(1,length(node_x));
R = zeros(1,sat_num);
delta_phi = zeros(1,sat_num);
signal_receive = zeros(1,length(t));
%%
   for m = 1:length(node_x)
        for k = 1:sat_num
            R(k) = sqrt((y(node_y(m))-sat_loc_y(k)).^2+(x(node_x(m))-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            signal_receive = signal_receive + s(k,:)*exp(1i*delta_phi(k));
        end
        cor(m) = signal_receive*signal_receive'/length(signal_receive);
        signal_receive = zeros(1,length(t));
   end
% 创建坐标网格
% [X, Y] = meshgrid(x, y);
% 
% % 初始化距离和相位差矩阵
% R = zeros(length(x), length(y), sat_num);
% delta_phi = zeros(length(x), length(y), sat_num);
% 
% % 计算距离和相位差矩阵
% for k = 1:sat_num
%     R(:, :, k) = sqrt((Y - sat_loc_y(k)).^2 + (X - sat_loc_x(k)).^2 + D^2);
%     delta_phi(:, :, k) = R(:, :, k) * 2 * pi / lambda;
% end
% 
% % 计算接收信号矩阵
% signal_receive = zeros(length(x), length(y), length(t));
% for k = 1:sat_num
%     signal_receive = signal_receive + reshape(s(k, :), 1, 1, []) .* exp(1i * delta_phi(:, :, k));
% end
% 
% % 计算相关矩阵
% cor = sum(signal_receive .* conj(signal_receive), 3) / length(t);

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
        score = score + cor(m);
end
score = -score/length(node_x);
end