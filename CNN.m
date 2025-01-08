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
idx = size(node,2);
%%
% node = node(:,1:round(idx/10));
% parameters = parameters(:,1:round(idx/10));
% 
% idx = size(node,2);

par = 1.25;%all/par = 训练集比例
train_x = node(:,1:round(idx/par)); 
train_y = parameters(:,1:round(idx/par));
test_x = node(:,round(idx/par):idx);
test_y = parameters(:,round(idx/par):idx);

method=@mapminmax;
%method=@mapstd;
[train_x,input_ps]=method(train_x);
test_x=method('apply',test_x,input_ps);
[train_y,output_ps]=method(train_y);
test_y=method('apply',test_y,output_ps);

train_x = reshape(train_x,10,2,size(train_x,2));
test_x = reshape(test_x,10,2,size(test_x,2));
idx_train = size(train_x,3);
idx_test = size(test_x,3);
degree_x = size(test_x,1);
degree_y = size(test_y,1);

input_train=reshape(train_x,[degree_x,2,1,idx_train]);%训练集输入
input_test=reshape(test_x,[degree_x,2,1,idx_test]);%测试集输入
output_train = train_y;%训练集输出
output_test  = test_y;%测试集输出
%%
load('sat4_phiAloc_model.mat');
close all
layers = [
    imageInputLayer([10 2 1])
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(8)
    regressionLayer];
    %'LearnRateSchedule',"piecewise"
options = trainingOptions('adam', ...
    'MaxEpochs',20,...
    'MiniBatchSize',16,...
    'InitialLearnRate',0.001,... 
    'GradientThreshold',1, ...
    'Verbose',false,...
    'Plots','training-progress',...
    'ValidationData',{input_test,output_test'});
net = trainNetwork(input_train, output_train', layers, options);
save('sat4_phiAloc_model.mat', 'net', 'options', 'layers');

%%
tic
YPred = predict(net,input_test);
toc

predict_value=method('reverse',YPred',output_ps);
predict_value=double(predict_value);
true_value=method('reverse',output_test,output_ps);

sat_num = 4;
predict_phi = predict_value(1:sat_num,:);
predict_loc = predict_value(sat_num+1:sat_num*2,:);


true_phi = true_value(1:sat_num,:);
true_loc = true_value(sat_num+1:sat_num*2,:);


figure
plot(1:idx_test,true_phi(1,:),'-*','linewidth',3)
hold on
plot(1:idx_test,predict_phi(1,:),'-s','linewidth',3)
legend('实际相位','预测相位')
grid on

figure
plot(1:idx_test,true_loc(1,:),'-*','linewidth',3)
hold on
plot(1:idx_test,predict_loc(1,:),'-s','linewidth',3)
legend('实际位置','预测位置')
grid on


rmse=sqrt(mean((true_value-predict_value).^2));
disp(['根均方差(RMSE)：',num2str(rmse)])
mae=mean(abs(true_value-predict_value));
disp(['平均绝对误差（MAE）：',num2str(mae)])
mape=mean(abs((true_value-predict_value)./true_value));
disp(['平均相对百分误差（MAPE）：',num2str(mape*100),'%'])

%%
test = zeros(10,2);
test(:,1) = [117	70	147	151	184	156	110	166	173	170];
test(:,2) = [163	35	59	104	201	184	115	9	113	36];
testss=mapminmax(test);
x = predict(net,testss);
predictss=method('reverse',x',output_ps);