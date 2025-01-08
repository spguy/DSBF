clc
clear
close all
tic
% basic parameters
sat_num = 4;
H = 0.1;%地面范围
step = H/201;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;
node_num = 10;
%node_x = randi(length(x_ground),1,node_num);
%node_y = randi(length(y_ground),1,node_num);
node_x = [117	70	147	151	184	156	110	166	173	170];
node_y = [ 163	35	59	104	201	184	115	9	113	36];
d = 1.92;%satellite distance
D = 550;%distance to the ground
sat_loc_x = [d d -d -d];
sat_loc_y = [d -d d -d];
nvars = sat_num*2;%变量个数

% 固定卫星功率和 
baselinex = [zeros(1,sat_num) ones(1,sat_num)];
[threshold,] = fitness_power(baselinex,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D);%以无优化的平均功率作为门限

%测试无优化时统计意义下的平均功率和门限作为参考
%[power_test,propability_test] = node_test(baselinex,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_num,D,threshold)
%%
baseline_power_rec = [];
best_power_rec = [];
baseline_propability_rec = [];
best_propability_rec = [];
for kk = -1:-0.5:-16%统计不同门限下的表现
threshold = kk;
disp('门限为')
disp(kk)

% 问题设置
%disp("以平均功率作为优化目标")
%targetFunction = @(x)fitness_power_nodeonly(x,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D); % 目标函数：调用 fitness 函数
%disp("以概率作为优化目标")
%targetFunction = @(x)fitness_propability_nodeonly(x,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D,threshold); % 目标函数：调用 fitness 函数
%disp("以功率和概率综合作为优化目标")
%targetFunction = @(x)fitness_multiobj_nodeonly(x,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D,threshold); % 目标函数：调用 fitness 函数
%disp("以加权综合作为优化目标")
targetFunction = @(x)fitness_weighted_nodeonly(x,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D,threshold); % 目标函数：调用 fitness 函数

constraintFunction = @(x)constraint(x,sat_num);

UB =[pi/2*ones(1,sat_num) 2*ones(1,sat_num)];
LB =[-pi/2*ones(1,sat_num) zeros(1,sat_num)];


% 使用 ga 函数进行单目标优化
%options = optimoptions('ga','Display','iter');
[x, fval, exitflag, output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,constraintFunction);

% 使用 gamultiobj 函数进行多目标优化
%options = optimoptions('gamultiobj','Display','iter');
%[x, fval, exitflag, output, population, scores] = gamultiobj(targetFunction, nvars, [], [], [], [], LB, UB,constraintFunction,options);

% 提取最优解中的 phi 和 A
optimalPhi = x(1:sat_num)/pi;
optimalA = x((sat_num+1):2*sat_num);

% 打印最终结果
disp(['最优解的相位值(pi)：', num2str(optimalPhi)])%pi
disp(['最优解的幅度：', num2str(optimalA)]);
disp(['约束值：', num2str(constraintFunction(x))]);

%% figures
[best_power,cor] = fitness_power(x(1,:),sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D);
best_propability = fitness_propability_nodeonly(x(1,:),sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D,threshold);
[baseline_power,cor_baseline] = fitness_power(baselinex,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D);
baseline_propability = fitness_propability_nodeonly(baselinex,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_x,node_y,D,threshold);
powerdistribution(cor,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
%powerdistribution(cor_baseline,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
sprintf("最优平均功率%f 无优化平均功率%f",best_power,baseline_power)
sprintf("最优概率%f 无优化概率%f",best_propability,baseline_propability)
baseline_power_rec = [baseline_power_rec,baseline_power];
best_power_rec = [best_power_rec,best_power];
baseline_propability_rec = [baseline_propability_rec,baseline_propability];
best_propability_rec = [best_propability_rec,best_propability];
end
%% pareto figures
% figure
% title("pareto fronts and baseline")
% scatter(fval(:,1),fval(:,2),'filled', 'b');
% hold on
% scatter(baseline_power,baseline_propability,'r', 'Marker', 'o');
% xlabel('平均功率')
% ylabel('概率')
% grid on
% toc
%%
figure
threshold_plot = 1:0.5:16;
plot(threshold_plot,-best_propability_rec,threshold_plot,-baseline_propability_rec,'LineWidth',3);
legend('优化后概率','优化前概率','Fontsize',16)
title('随门限增大的概率曲线','Fontsize',16)
grid on
figure
plot(threshold_plot,-best_power_rec,threshold_plot,-baseline_power_rec,'LineWidth',3)
legend('优化后平均功率','优化前平均功率','Fontsize',16)
title('随门限增大的平均功率曲线','Fontsize',16)
grid on