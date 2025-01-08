clc
clear
close all
tic
% basic parameters
sat_num = 4;

H = 1;%地面范围
step = H/201;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;

H_sat = 5;
step_sat = H_sat/201;
sat_available = 10;
y_sat = -H_sat/2:step_sat:H_sat/2;
x_sat = -H_sat/2:step_sat:H_sat/2;
sattochoose_x = x_sat(randperm(length(x_sat), sat_available));
sattochoose_y = y_sat(randperm(length(y_sat), sat_available));

node_num = 10;
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
nvars = sat_num*2+sat_num;%变量个数

% 固定卫星功率和 
basesat =  randi([1, sat_available], 1, 4);
baselinex = [zeros(1,sat_num) ones(1,sat_num) basesat];%[phi A sat_X sat_Y]
[baseline_power] = fitness_power_nodeonly_phiAloc_select(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,sattochoose_x,sattochoose_y);%以无优化的平均功率作为门限

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
targetFunction = @(x)fitness_weighted_nodeonly_phiAloc_select(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sattochoose_x,sattochoose_y); % 目标函数：调用 fitness 函数

constraintFunction = @(x)constraint(x,sat_num);

UB =[pi*ones(1,sat_num) ones(1,sat_num) sat_available/4 sat_available/2 sat_available*3/4 sat_available];%卫星分别位于四个象限
LB =[-pi*ones(1,sat_num) zeros(1,sat_num) 1 sat_available/4+1 sat_available/2+1 sat_available*3/4+1];


% 使用 ga 函数进行单目标优化
%options = optimoptions('ga','Display','iter');
%[x, fval, exitflag, output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,constraintFunction,options);
%[x, fval, exitflag, output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,[],options);

%使用 pso 函数进行单目标优化
options = optimoptions('particleswarm','Display','iter');
[x, fval, exitflag, output] = particleswarm(targetFunction, nvars, LB, UB,options);
toc
% 提取最优解中的 phi 和 A
optimalPhi = x(1:sat_num)/pi;
optimalA = x((sat_num+1):2*sat_num);
optimalidx = round(x(sat_num*2+1:sat_num*3));  


optimalXloc = sattochoose_x(optimalidx);
optimalYloc = sattochoose_y(optimalidx);
% 打印最终结果
disp(['最优解的相位值(pi)：', num2str(optimalPhi)])%pi
disp(['最优解的幅度：', num2str(optimalA)]);
sat_locations = [];
for i =1:sat_num
sat_locations = [sat_locations optimalXloc(i) optimalYloc(i)];
end
disp(['最优解的卫星位置：',num2str(sat_locations)])

%% figures
[t,best_power,best_propability] = fitness_weighted_nodeonly_phiAloc_select(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sattochoose_x,sattochoose_y);
[sxa,cor] = fitness_power_phiAloc_select(x,sat_num,x_ground,y_ground,node_x,node_y,D,sattochoose_x,sattochoose_y); 
[sxa,cor_baseline] = fitness_power_phiAloc_select(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,sattochoose_x,sattochoose_y); 
[a,b,baseline_propability] = fitness_weighted_nodeonly_phiAloc_select(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sattochoose_x,sattochoose_y);
%powerdistribution(cor,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
%powerdistribution(cor_baseline,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
sprintf("最优平均功率%f 无优化平均功率%f",best_power,baseline_power)
sprintf("最优概率%f 无优化概率%f",best_propability,baseline_propability)
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
scatter(optimalXloc, optimalYloc,'b');
hold on
scatter(sattochoose_x(basesat),sattochoose_y(basesat),'r')
legend('所有可选卫星','随机选取','通过算法选取')
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
A = xx((sat_num+1):sat_num*2);
satidx = xx(sat_num*2+1:sat_num*3);     

sat_loc_x = x_sat(round(satidx));
sat_loc_y = y_sat(round(satidx));

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











 
function [score,power,propability] = fitness_weighted_nodeonly_phiAloc_select(xx,sat_num,x,y,node_x,node_y,D,threshold,x_sat,y_sat)

%disp(['Debug: x = ', num2str(x)]);
phi = xx(1:sat_num);
A = xx((sat_num+1):sat_num*2);
satidx = xx(sat_num*2+1:sat_num*3);     

sat_loc_x = x_sat(round(satidx));
sat_loc_y = y_sat(round(satidx));

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


function [score,cor] = fitness_power_phiAloc_select(xx,sat_num,x,y,node_x,node_y,D,x_sat,y_sat)
%disp(['Debug: x = ', num2str(x)]);
phi = xx(1:sat_num);
A = xx((sat_num+1):sat_num*2);
satidx = xx(sat_num*2+1:sat_num*3);     

sat_loc_x = x_sat(round(satidx));
sat_loc_y = y_sat(round(satidx));

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