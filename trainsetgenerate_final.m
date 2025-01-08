clc 
clear
close all
% basic parameters
sat_num = 4;
node_num = 10;
node_loc = zeros(node_num*2,1);
philoc = zeros(sat_num*2,1);
H = 1;%地面范围
step = H/201;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;    


H_sat = 10;
step_sat = H_sat/201;
sat_available = 10;
y_sat = -H_sat/2:step_sat:H_sat/2;
x_sat = -H_sat/2:step_sat:H_sat/2;
sattochoose_x = x_sat(randperm(length(x_sat), sat_available));%可选卫星的坐标
sattochoose_y = y_sat(randperm(length(y_sat), sat_available));

D = 550;%distance to the ground
nvars = sat_num*2;%变量个数
d = 1.92; 

tic
for i = 1:5000
node_x = randi(length(x_ground),1,node_num);
node_y = randi(length(y_ground),1,node_num);
nodes = [x_ground(node_x),y_ground(node_y)]';
% 固定卫星功率
basesat =  randi([1, sat_available], 1, 4);%随机选星作为baseline
baselinex = [zeros(1,sat_num) basesat];
[baseline_power] = fitness_power_nodeonly_phiAloc_select(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,sattochoose_x,sattochoose_y);%以无优化的平均功率作为门限
threshold = baseline_power;
targetFunction = @(x)fitness_weighted_nodeonly_phiAloc_select(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold,sattochoose_x,sattochoose_y); % 目标函数：调用 fitness 函数

UB =[pi*ones(1,sat_num) sat_available/4 sat_available/2 sat_available*3/4 sat_available];%为了避免重复，将选星分为四个区间，每个变量从一个区间中选取
LB =[-pi*ones(1,sat_num) 1 sat_available/4 sat_available/2 sat_available*3/4];

%options = optimoptions('particleswarm','Display','iter');
[x, fval, exitflag, output] = particleswarm(targetFunction, nvars, LB, UB);
node_loc = cat(2,node_loc,nodes);
philoc = cat(2,philoc,x');
end
toc

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

function [score,power,propability] = fitness_weighted_nodeonly_phiAloc_select(xx,sat_num,x,y,node_x,node_y,D,threshold,x_sat,y_sat)

%disp(['Debug: x = ', num2str(x)]);
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
