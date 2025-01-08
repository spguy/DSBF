clc
clear
close all

% basic parameters
sat_num = 4;
H = 1;%地面范围
step = H/201;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;
node_num = 10;
node_x = randi(length(x_ground),1,node_num);
node_y = randi(length(y_ground),1,node_num);
%好靶子
%node_x = [1907	1778	88	5	1429	1149	504	697	1593	1908];
%node_y = [ 1501	391	1601	1999	102	650	721	976	1838	378];
d = 19.2;%satellite distance
D = 550;%distance to the ground
nvars = sat_num*2+sat_num;%变量个数
%% 产生四条轨道
sat_track_theta = 0 + rand(1,sat_num) * (180 - 0);%0到180度内的随机数，四条轨道
sat_track_bias = -d + rand(1,sat_num) * (d - (-d));%-d到d之间的截距
%sat_track_slope = [];
% 固定卫星功率和 
baselinex = [zeros(1,sat_num) ones(1,sat_num) d d d d];%[phi A 四个卫星的横坐标]
[baseline_power] = fitness_power_nodeonly_phiAloc_track(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,sat_track_theta,sat_track_bias);%以无优化的平均功率作为门限
threshold = -3;
%测试无优化时统计意义下的平均功率和门限作为参考
%[power_test,propability_test] = node_test(baselinex,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_num,D,threshold)
%%
baseline_power_rec = [];
best_power_rec = [];
baseline_propability_rec = [];
best_propability_rec = [];
% for kk = 0:-0.5:-sat_num.^2
%     threshold = kk;

% 问题设置
targetFunction = @(x)fitness_weighted_nodeonly_phiAloc_track(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sat_track_theta,sat_track_bias); % 目标函数：调用 fitness 函数

constraintFunction = @(x)constraint(x,sat_num);

UB =[pi*ones(1,sat_num) ones(1,sat_num) d*ones(1,sat_num)];%卫星最大半径为d
LB =[-pi*ones(1,sat_num) zeros(1,sat_num) -d*ones(1,sat_num)];

tic;
% 使用 ga 函数进行单目标优化
options = optimoptions('particleswarm');

[x, fval, exitflag, output] = particleswarm(targetFunction, nvars, LB, UB,options);
calculation_time = toc;
% 提取最优解中的 phi 和 A
optimalPhi = x(1:sat_num)/pi;
optimalA = x((sat_num+1):2*sat_num);
optimalXloc = x(sat_num*2+1:sat_num*3);  
optimalYloc = x(sat_num*2+1:sat_num*3).*tand(sat_track_theta)+sat_track_bias;

% 打印最终结果
disp(['最优解的相位值(pi)：', num2str(optimalPhi)])%pi
disp(['最优解的幅度：', num2str(optimalA)]);
sat_locations = [];
for i =1:sat_num
sat_locations = [sat_locations optimalXloc(i) optimalYloc(i)];
end
disp(['最优解的卫星位置：',num2str(sat_locations)])
disp(['约束值：', num2str(constraintFunction(x))]);

%% figures
[t,best_power,best_propability] = fitness_weighted_nodeonly_phiAloc_track(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sat_track_theta,sat_track_bias);
%%
sat_speed = 7.8; %卫星速度为7.8km/s = 第一宇宙速度
sat_moving_x = sat_speed * calculation_time;%卫星在计算过程中的移动距离
x_moving = x;
x_moving(sat_num*2+1:sat_num*3) = x_moving(sat_num*2+1:sat_num*3) + sat_moving_x.*cosd(sat_track_theta);
node_x_random = zeros(1,length(node_x));
node_y_random = zeros(1,length(node_y));
for a = 1:length(node_x)
    node_random = randi([-1,1]);
    node_x_random(a) = node_x(a) + node_random;
    node_y_random(a) = node_y(a) + node_random;
end

[tt,best_power_moving,best_propability_moving] = fitness_weighted_nodeonly_phiAloc_track(x_moving,sat_num,x_ground,y_ground,node_x_random,node_y_random,D,threshold,sat_track_theta,sat_track_bias);

[sxa,cor] = fitness_power_phiAloc_track(x,sat_num,x_ground,y_ground,node_x,node_y,D,sat_track_theta,sat_track_bias); 
[sxa,cor_baseline] = fitness_power_phiAloc_track(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,sat_track_theta,sat_track_bias); 
[a,b,baseline_propability] = fitness_weighted_nodeonly_phiAloc_track(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sat_track_theta,sat_track_bias);
powerdistribution(cor,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
powerdistribution(cor_baseline,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
sprintf("最优平均功率%f 无优化平均功率%f",best_power,baseline_power)
sprintf("最优概率%f 无优化概率%f",best_propability,baseline_propability)
sprintf("考虑卫星位移和地面端估计不准后的优化功率%f 和 概率%f",best_power_moving,best_propability_moving)
baseline_power_rec = [baseline_power_rec,baseline_power];
best_power_rec = [best_power_rec,best_power];
baseline_propability_rec = [baseline_propability_rec,baseline_propability];
best_propability_rec = [best_propability_rec,best_propability];
disp(calculation_time)
%end
%%
% figure
% threshold_plot = 0:0.5:16;
% plot(threshold_plot,-best_propability_rec,threshold_plot,-baseline_propability_rec,'LineWidth',3);
% legend('优化后概率','优化前概率','Fontsize',16)
% title('随门限增大的概率曲线','Fontsize',16)
% grid on
% figure
% plot(threshold_plot,-best_power_rec,threshold_plot,-baseline_power_rec,'LineWidth',3)
% legend('优化后平均功率','优化前平均功率','Fontsize',16)
% title('随门限增大的平均功率曲线','Fontsize',16)
% grid on
% 
%%
figure
scatter(optimalXloc,optimalYloc)
% 添加坐标轴
xlim([-2*d,2*d]);
ylim([-2*d,2*d]);
hold on;
plot([0, 0], ylim, 'k--'); % 垂直坐标轴
plot(xlim, [0, 0], 'k--'); % 水平坐标轴
title('卫星轨道和卫星坐标')
grid on
% 把所有点标红
scatter(optimalXloc, optimalYloc,50, 'b','filled');
for i =1:sat_num
plot(xlim,tand(sat_track_theta(i)).*xlim+sat_track_bias(i))
hold on
end
%%
function [score] = fitness_power_nodeonly_phiAloc_track(xx,sat_num,x,y,node_x,node_y,D,sat_track_theta,sat_track_bias)
phi = xx(1:sat_num);
A = xx((sat_num+1):sat_num*2);
sat_loc_x = xx(sat_num*2+1:sat_num*3);     
sat_loc_y = xx(sat_num*2+1:sat_num*3).*tand(sat_track_theta)+sat_track_bias;
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
score = 0;
for m = 1:length(node_x)
        score = score + cor(m);
end
score = -score/length(node_x);
end











 
function [score,power,propability] = fitness_weighted_nodeonly_phiAloc_track(xx,sat_num,x,y,node_x,node_y,D,threshold,sat_track_theta,sat_track_bias)

%disp(['Debug: x = ', num2str(x)]);
phi = xx(1:sat_num);
A = xx((sat_num+1):sat_num*2);
sat_loc_x = xx(sat_num*2+1:sat_num*3);     
sat_loc_y = xx(sat_num*2+1:sat_num*3).*tand(sat_track_theta)+sat_track_bias;
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

%threshold = 4; 
score_power = 0;
score_propability = 0;
for m = 1:length(node_x)
    score_power = score_power + cor(m);
    if cor(m)>abs(threshold)
        score_propability = score_propability + 1;
    end
end
score_power = -score_power/length(node_x);%求平均
power = score_power;
score_power = score_power/sat_num.^2;%归一化
score_propability = -score_propability/length(node_x);
propability = score_propability;
weight_power = 0.5;
weight_propability = 1-weight_power;
score = weight_power*score_power + weight_propability*score_propability;
end


function [score,cor] = fitness_power_phiAloc_track(xx,sat_num,x,y,node_x,node_y,D,sat_track_theta,sat_track_bias)
%disp(['Debug: x = ', num2str(x)]);
phi = xx(1:sat_num);
A = xx((sat_num+1):sat_num*2);
sat_loc_x = xx(sat_num*2+1:sat_num*3);     
sat_loc_y = xx(sat_num*2+1:sat_num*3).*tand(sat_track_theta)+sat_track_bias;
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
            signal_receive = signal_receive + s(k,:)*exp(1i*delta_phi(k));
        end
        cor(m,n) = signal_receive*signal_receive'/length(signal_receive);
        signal_receive = zeros(1,length(t));
    end
end

score = 0;
for m = 1:length(node_x)
        score = score + cor(node_x(m),node_y(m));
end
score = -score/length(node_x);
end