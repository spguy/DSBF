clc
clear
close all
% 优化的变量：选星（从给定的sat_available = 10个卫星中选取sat_num = 4个卫星）； 每颗卫星的相位phi
% 解的构成：[phi*sat_num index*sat_num(选取的卫星标号)]
% 目标函数：地面端节点的平均功率和门限概率的权重符合函数
% 优化算法：ga/pso
% 选星区域：直径H_sat = 10km的范围内 
% **为了避免选取重复的卫星，思路1：将所有可选卫星分为4个区域，选取在区域内发生 
% **避免重复思路2：当选取到重复卫星后，检索物理距离最近的卫星作为替代
% 地面端用户分布范围： 直径H = 1km内的node_num = 10个用户
% 轨道高度 D = 550km

% 计算选星和相位
% 根据轨道计算卫星位置和相位
% 补偿该相位变化
% basic parameters

calculation_time = 0;%初始化计算时间

sat_num = 4;

H = 1;%地面范围
step = H/2001;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;

H_sat = 10; %卫星分布范围
step_sat = H_sat/201;
sat_available = 20;
v_sat = 7.8;%卫星运动速度为7.8km/s
y_sat = -H_sat/2:step_sat:H_sat/2;
x_sat = -H_sat/2:step_sat:H_sat/2;
% 产生可选卫星的坐标
sattochoose_x = x_sat(randperm(length(x_sat), sat_available));%可选卫星的坐标
%sattochoose_x = [-0.422885572139303	-3.85572139303483	-2.96019900497512	1.66666666666667	3.35820895522388	-3.35820895522388	-2.81094527363184	-3.60696517412935	3.15920398009950	-4.20398009950249	3.40796019900498	-1.46766169154229	-4.95024875621891	-4.25373134328358	2.21393034825871	0.621890547263681	-1.71641791044776	-0.174129353233830	-1.51741293532338	0.273631840796019];
sattochoose_y = y_sat(randperm(length(y_sat), sat_available));
%sattochoose_y = [2.26368159203980	-4.05472636815920	-1.66666666666667	1.06965174129353	-4.70149253731343	-3.55721393034826	2.46268656716418	0.124378109452736	4.35323383084577	3.20895522388060	-0.970149253731343	-0.323383084577114	-2.51243781094527	1.71641791044776	0.273631840796019	2.91044776119403	3.90547263681592	0.920398009950248	1.41791044776119	1.01990049751244];
% 产生可选卫星的轨道倾角（初始位置为随机生成的位置），用来计算卫星运动轨迹
sattrack_theta =  0 + rand(1,sat_available) * (180 - 0);%0到180度内的随机数

node_num = 3;
node_x = randi(length(x_ground),1,node_num);
node_y = randi(length(y_ground),1,node_num);
%好靶子
%node_x = [1907	1778	88	5	1429	1149	504	697	1593	1908];

%node_y = [ 1501	391	1601	1999	102	650	721	976	1838	378];
%对比组
%node_x = [117	70	147	151	184	156	110	166	173	170];
%node_y = [ 163	35	59	104	201	184	115	9	113	36];
d = 1.92;%satellite distance
D = 550;%distance to the ground
nvars = sat_num*2;%变量个数

% 固定卫星功率和 
%basesat =  randi([1, sat_available], 1, 4);%随机选星作为baseline
basesat = [13 11 20 12];
baselinex = [zeros(1,sat_num) basesat];
[baseline_power] = fitness_power_nodeonly_phiAloc_select(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,sattochoose_x,sattochoose_y);%以无优化的平均功率作为门限

threshold = baseline_power;
%测试无优化时统计意义下的平均功率和门限作为参考
%[power_test,propability_test] = node_test(baselinex,sat_num,sat_loc_x,sat_loc_y,x_ground,y_ground,node_num,D,threshold)
%%
baseline_power_rec = [];
best_power_rec = [];
baseline_propability_rec = [];
best_propability_rec = [];
% for kk = 0:-0.5:-sat_num.^2 用来统计门限变化的目标函数值曲线
%     threshold = kk;

% 问题设置
moving0 = 0;%以当前卫星位置进行计算
moving1 = 1;%考虑卫星运动
targetFunction = @(x)fitness_weighted_nodeonly_phiAloc_select(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sattochoose_x,sattochoose_y,sattrack_theta,v_sat,calculation_time,moving0); % 目标函数：调用 fitness 函数

%constraintFunction = @(x)constraint(x,sat_num);

UB =[pi*ones(1,sat_num) sat_available/4 sat_available/2 sat_available*3/4 sat_available];%为了避免重复，将选星分为四个区间，每个变量从一个区间中选取
LB =[-pi*ones(1,sat_num) 1 sat_available/4 sat_available/2 sat_available*3/4];

%UB =[pi*ones(1,sat_num) sat_available sat_available sat_available sat_available];%不加限制，后续处理避免重复
%LB =[-pi*ones(1,sat_num) 1 1 1 1];

tic;
% 使用 ga 函数进行单目标优化
options = optimoptions('ga','PopulationSize',50,'MaxGenerations',10);%default: ps 200; mg 200*number of variable
%[x, fval, exitflag, output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,constraintFunction,options);
[x, fval, exitflag, output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,[],options);

%使用 pso 函数进行单目标优化
%options = optimoptions('particleswarm','Display','iter');
%[x, fval, exitflag, output] = particleswarm(targetFunction, nvars, LB, UB,options);

calculation_time = toc

% 提取最优解中的 phi 和 loc
optimalPhi = x(1:sat_num)/pi;
optimalidx = round(x(sat_num+1:sat_num*2));  
optimalXloc = sattochoose_x(optimalidx);
optimalYloc = sattochoose_y(optimalidx);
% 打印最终结果
disp(['最优解的相位值(pi)：', num2str(optimalPhi)])%pi
sat_locations = [];
for i =1:sat_num
sat_locations = [sat_locations optimalXloc(i) optimalYloc(i)];
end
disp(['最优解的卫星位置：',num2str(sat_locations)])

%% 计算优化结果（需要考虑卫星的运动）

[t,best_power,best_propability,sat_loc_x_moving,sat_loc_y_moving] = fitness_weighted_nodeonly_phiAloc_select(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sattochoose_x,sattochoose_y,sattrack_theta,v_sat,calculation_time,moving1);
[t,best_power_m0,best_propability_m0] = fitness_weighted_nodeonly_phiAloc_select(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sattochoose_x,sattochoose_y,sattrack_theta,v_sat,calculation_time,moving0);
[sxa,cor] = fitness_power_phiAloc_select(x,sat_num,x_ground,y_ground,node_x,node_y,D,sattochoose_x,sattochoose_y,sattrack_theta,v_sat,calculation_time,moving1); 
[sxa,cor_baseline] = fitness_power_phiAloc_select(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,sattochoose_x,sattochoose_y,sattrack_theta,v_sat,calculation_time,moving0); 

[a,b,baseline_propability] = fitness_weighted_nodeonly_phiAloc_select(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sattochoose_x,sattochoose_y,sattrack_theta,v_sat,calculation_time,moving0);
%powerdistribution(cor,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
%powerdistribution(cor_baseline,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
sprintf("最优平均功率%f 无优化平均功率%f",best_power,baseline_power)
sprintf("最优概率%f 无优化概率%f",best_propability,baseline_propability)
sprintf("不考虑移动的最优平均功率%f 最优概率%f",best_power_m0,best_propability_m0)
% baseline_power_rec = [baseline_power_rec,baseline_power];
% best_power_rec = [best_power_rec,best_power];
% baseline_propability_rec = [baseline_propability_rec,baseline_propability];
% best_propability_rec = [best_propability_rec,best_propability];
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
figure
scatter(sattochoose_x,sattochoose_y,'x')
hold on
scatter(optimalXloc, optimalYloc,50,'b');
hold on
scatter(sattochoose_x(basesat),sattochoose_y(basesat),'r')
legend('所有可选卫星','通过算法选取','随机选取')
% 添加坐标轴
% xlim([-H_sat/2,H_sat/2]);
% ylim([-H_sat/2,H_sat/2]);
% hold on;
% plot([0, 0], ylim, 'k--'); % 垂直坐标轴
% plot(xlim, [0, 0], 'k--'); % 水平坐标轴
grid on
% 把所有点标红


%%
function [score] = fitness_power_nodeonly_phiAloc_select(xx,sat_num,x,y,node_x,node_y,D,x_sat,y_sat)
phi = xx(1:sat_num);
satidx = xx(sat_num+1:sat_num*2);     
sat_loc_x = x_sat(round(satidx));
sat_loc_y = y_sat(round(satidx));
A = ones(1,sat_num);
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
%satellite coordinates 
t = 0:1e-10:1e-9;
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











 
function [score,power,propability,x_moving,y_moving] = fitness_weighted_nodeonly_phiAloc_select(xx,sat_num,x,y,node_x,node_y,D,threshold,x_sat,y_sat,sattrack_theta,v_sat,calculation_time,moving)

%disp(['Debug: x = ', num2str(x)]);
phi = xx(1:sat_num);
satidx = xx(sat_num+1:sat_num*2);
% 计算运算结束时卫星的位置
x_moving = moving.*v_sat.*calculation_time.*cosd(sattrack_theta(round(satidx)));
y_moving = moving.*v_sat.*calculation_time.*sind(sattrack_theta(round(satidx)));
sat_loc_x = x_sat(round(satidx))+x_moving;
sat_loc_y = y_sat(round(satidx))+y_moving;

A = ones(1,sat_num);
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
%satellite coordinates 
t = 0:1e-10:1e-9;
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


function [score,cor] = fitness_power_phiAloc_select(xx,sat_num,x,y,node_x,node_y,D,x_sat,y_sat,sattrack_theta,v_sat,calculation_time,moving)
%disp(['Debug: x = ', num2str(x)]);
phi = xx(1:sat_num);
satidx = xx(sat_num+1:sat_num*2);     
sat_loc_x = x_sat(round(satidx))+moving*v_sat*calculation_time*cosd(sattrack_theta(round(satidx)));
sat_loc_y = y_sat(round(satidx))+moving*v_sat*calculation_time*sind(sattrack_theta(round(satidx)));
A = ones(1,sat_num);
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
%satellite coordinates 
t = 0:1e-10:1e-9;
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