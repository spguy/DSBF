clc
clear
close all
tic
% basic parameters
sat_num = 4;
H = 1;%地面范围
step = H/2001;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;
node_num = 10;
%node_x = randi(length(x_ground),1,node_num);
%node_y = randi(length(y_ground),1,node_num);
%好靶子
node_x = [1907	1778	88	5	1429	1149	504	697	1593	1908];
node_y = [ 1501	391	1601	1999	102	650	721	976	1838	378];
%对比组
%node_x = [117	70	147	151	184	156	110	166	173	170];
%node_y = [ 163	35	59	104	201	184	115	9	113	36];
d = 1.92;%satellite distance
D = 550;%distance to the ground
nvars = sat_num*2+sat_num*2;%变量个数

% 固定卫星功率和 
baselinex = [zeros(1,sat_num) ones(1,sat_num) d d -d -d d -d d -d];%[phi A sat_X sat_Y]
[baseline_power] = fitness_power_nodeonly_phiAloc(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D);%以无优化的平均功率作为门限
threshold = baseline_power;
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
targetFunction = @(x)fitness_weighted_nodeonly_phiAloc(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold); % 目标函数：调用 fitness 函数

constraintFunction = @(x)constraint(x,sat_num);

UB =[pi*ones(1,sat_num) ones(1,sat_num) d 0 0 d d d 0 0];%卫星分别位于四个象限
LB =[-pi*ones(1,sat_num) zeros(1,sat_num) 0 -d -d 0 0 0 -d -d];


% 使用 pso 函数进行单目标优化
options = optimoptions('particleswarm','Display','iter');

[x, fval, exitflag, output] = particleswarm(targetFunction, nvars, LB, UB,options);
toc
% 提取最优解中的 phi 和 A
optimalPhi = x(1:sat_num)/pi;
optimalA = x((sat_num+1):2*sat_num);
optimalXloc = x(sat_num*2+1:sat_num*3);  
optimalYloc = x(sat_num*3+1:sat_num*4);

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
[t,best_power,best_propability] = fitness_weighted_nodeonly_phiAloc(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold);
[sxa,cor] = fitness_power_phiAloc(x,sat_num,x_ground,y_ground,node_x,node_y,D); 
[sxa,cor_baseline] = fitness_power_phiAloc(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D); 
[a,b,baseline_propability] = fitness_weighted_nodeonly_phiAloc(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,threshold);
powerdistribution(cor,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
powerdistribution(cor_baseline,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
sprintf("最优平均功率%f 无优化平均功率%f",best_power,baseline_power)
sprintf("最优概率%f 无优化概率%f",best_propability,baseline_propability)
baseline_power_rec = [baseline_power_rec,baseline_power];
best_power_rec = [best_power_rec,best_power];
baseline_propability_rec = [baseline_propability_rec,baseline_propability];
best_propability_rec = [best_propability_rec,best_propability];
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
% 
% figure
% scatter(optimalXloc,optimalYloc)
% % 添加坐标轴
% xlim([-d,d]);
% ylim([-d,d]);
% hold on;
% plot([0, 0], ylim, 'k--'); % 垂直坐标轴
% plot(xlim, [0, 0], 'k--'); % 水平坐标轴
% grid on
% % 把所有点标红
% scatter(optimalXloc, optimalYloc,20, 'b','filled');