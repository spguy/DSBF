%%%%%%%%%%%%%%%%%%%%%输出的编码为[phi A sat_X sat_Y]%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
warning off;
node_num = 10;
addpath(genpath(pwd));
rng('default')
node_loc = load('C:\Users\pc\Desktop\训练集\nodes.mat');
philoc = load('C:\Users\pc\Desktop\训练集\philoc.mat');
node = node_loc.node_loc;
parameters = philoc.philoc;
%%
% idx = size(node,2);
% node = node(:,1:round(idx/10));
% parameters = parameters(:,1:round(idx/10));


method=@mapminmax;
%method=@mapstd;
[node input_ps]= method(node);
[parameters output_ps]=method(parameters);

%%
% 创建神经网络
net = feedforwardnet(50);  % 50个隐藏层神经元
net = configure(net,node,parameters);
% 训练神经网络
[trainednet,tr] = train(net,node,parameters);

%%
% sat_num = 4;
% H = 0.1;%地面范围
% step = H/201;
% y_ground = -H/2:step:H/2;
% x_ground = -H/2:step:H/2;
% node_num = 10;
% node_x = [117	70	147	151	184	156	110	166	173	170];
% node_y = [ 163	35	59	104	201	184	115	9	113	36];
% d = 1.92;%satellite distance
% D = 550;%distance to the ground
% nvars = sat_num*2+sat_num*2;%变量个数
% YPred = net([node_x node_y]');
% predict_value = method_y('reverse',YPred,output_ps);
% baselinex = [zeros(1,sat_num) ones(1,sat_num) d d -d -d d -d d -d];%[phi A sat_X sat_Y]
% [baseline_power] = fitness_power_nodeonly_phiAloc(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D);%以无优化的平均功率作为门限
% threshold = baseline_power;
% [t,best_power,best_propability] = fitness_weighted_nodeonly_phiAloc(predict_value,sat_num,x_ground,y_ground,node_x,node_y,D,threshold);
% [sxa,cor] = fitness_power_phiAloc(predict_value,sat_num,x_ground,y_ground,node_x,node_y,D); 
% [sxa,cor_baseline] = fitness_power_phiAloc(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D); 
% [a,b,baseline_propability] = fitness_weighted_nodeonly_phiAloc(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D,threshold);
% powerdistribution(cor,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
% powerdistribution(cor_baseline,x_ground,y_ground,node_x,node_y,x_ground,y_ground);
% sprintf("最优平均功率%f 无优化平均功率%f",best_power,baseline_power)
% sprintf("最优概率%f 无优化概率%f",best_propability,baseline_propability)
%%
YPred = net(node);


% predict_value = YPred;
% true_value = output_test;

% 反归一化
predict_value=method('reverse',YPred,output_ps);
predict_value=double(predict_value);
idx_test = size(parameters,2);

figure
plot(1:idx_test,parameters(1,:),'-*','linewidth',3)
hold on
plot(1:idx_test,predict_value(1,:),'-s','linewidth',3)
legend('实际相位','预测相位')
grid on
sat_num = 4;
figure
plot(1:idx_test,parameters(sat_num+1,:),'-*','linewidth',3)
hold on
plot(1:idx_test,predict_value(sat_num+1,:),'-s','linewidth',3)
legend('实际位置','预测位置')
grid on



