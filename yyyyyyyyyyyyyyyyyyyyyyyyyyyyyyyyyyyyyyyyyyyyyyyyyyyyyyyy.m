clc
clear
close all

%%%%%%%%%算法基本参数
% 优化的变量：选星（从给定的sat_available = 10个卫星中选取sat_num = 4个卫星）； 每颗卫星的相位phi
% 解的构成：[phi*sat_num index*sat_num(选取的卫星标号)]
% 目标函数：地面端单个节点的功率，考虑传输过程中的衰减
% 优化算法：ga/pso
% 选星区域：直径H_sat = 10km的范围内 
% 轨道高度 D = 550km

%%%%%%%%%算法基本流程
% 计算选星和相位
% 根据轨道计算卫星位置和相位
% 补偿该相位变化

% basic parameters
c = 3e8;
frequency = 3.5e9;
H_sat = 2*1000; %卫星分布范围
H = 1*1000;%地面范围
w = frequency*2*pi;
lambda = c/frequency;%波长
sat_num = 4;
v_sat = 7.8*1000;%卫星运动速度为7.8km/s
D = 600*1000;%orbit 
node_num = 7;%地面节点数
step = H/2001;%地面端粒度，只会影响地面端功率绘图的速度
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;
powerBaseExp = [];
powerGAExp = [];
proGAExp = [];
scoreGAExp = [];
for experiment = 1:20
%根据节点数随机产生地面节点
rng('default');
rng(experiment)%6 为较优情况

node_x = randi(length(x_ground),1,node_num);
node_y = randi(length(y_ground),1,node_num);
step_sat = H_sat/201;%卫星分布粒度，没什么太大影响
sat_available = 10;%可选卫星的数量
sat_x = -H_sat/2:step_sat:H_sat/2;
sat_y = -H_sat/2:step_sat:H_sat/2;
nvars = 2*sat_num;%优化目标为相位和四颗卫星的标号
%各类指标记录的初始化，必须放在循环开始前，用来记录整个流程的指标变化
score_record = [];%用来记录每次更新的得分
propability_record = [];%用来记录每次更新的地面概率
power_record = [];%用来记录每次更新的地面功率
%GA结果记录
score_recordGA = [];%用来记录每次更新的得分
propability_recordGA = [];%用来记录每次更新的地面概率
power_recordGA = [];%用来记录每次更新的地面功率
%baseline和持续时间记录
power_recordu0 = [];%用来记录权重u=0时的baseline
propability_recordu0 = [];
score_recordu0 = [];
outtime_record = [];%用来记录每次选星后的持续时间

calculationmin_record = [];
calculationGA_record = [];%用来对比计算速度
figure%卫星图，在循环开始前将图片打开

selectMax = 1;%决定重复选星几次

power_base = 8;%对标MIMO
%% 开始选星循环
%for a = 1:selectMax

%根据可选卫星数量随机产生可选卫星的坐标和每个卫星的轨道
sattochoose_x = sat_x(randperm(length(sat_x), sat_available));
sattochoose_y = sat_y(randperm(length(sat_y), sat_available));
sattrack_theta =  0 + rand(1,sat_available) * (90 - 0);%0到90度内的随机数

%随即选星，0相位作为基准
selectionNum = nchoosek(sat_available,sat_num);
selectionSatAll = [];
for i = 1:sat_available
    for j = i+1:sat_available
        for k = j+1:sat_available
            for p = k+1:sat_available
                selectionSatAll = [selectionSatAll;i j k p];
            end
        end
    end
end

x_base = [zeros(1,sat_num) randi(selectionNum,1,1)];
[fval_base] = score(x_base,sat_num,x_ground,y_ground,node_x,node_y,sattochoose_x,sattochoose_y,power_base,frequency,D,selectionSatAll);
powerBaseExp = [powerBaseExp -fval_base];
disp(['无优化分数：', num2str(-fval_base)])%pi

phii = x_base(1:sat_num);
satidxi = x_base(sat_num+1);  
satSelectioni = selectionSatAll(satidxi,:);
sat_loc_xi = sattochoose_x(satSelectioni);
sat_loc_yi = sattochoose_y(satSelectioni);
 



 
%% 画出随机选星，无相位情况下的地面分布
[aaa,corRandom] = power_calculation(sat_num,sat_loc_xi,sat_loc_yi,phii,x_ground,y_ground,H,node_x,node_y,D,frequency);
%powerdistribution(corRandom,x_ground,y_ground,node_x,node_y)
%%
%%%%%%%%%使用 ga 函数进行单目标优化

%输入变量长度为sat_num*2，对应相位和卫星标号
targetFunction = @(x)score(x,sat_num,x_ground,y_ground,node_x,node_y,sattochoose_x,sattochoose_y,power_base,frequency,D,selectionSatAll);

%为了避免重复，将选星分为四个区间，每个变量从一个区间中选取
UB =[pi*ones(1,sat_num)  selectionNum];
LB =[-pi*ones(1,sat_num)  1];

%设置options，减少populationsize和maxgeneration以节约计算时间
%options = optimoptions('ga','PopulationSize',50,'MaxGenerations',10);%default: ps 200; mg 200*number of variable
options = optimoptions('ga');

tic;
[x, fval,exitflag,output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,[],options);
calculation_time = toc%记录计算时间
[fval,power,propability] = score(x,sat_num,x_ground,y_ground,node_x,node_y,sattochoose_x,sattochoose_y,power_base,frequency,D,selectionSatAll);%算出最优情况下的各项指标
powerGAExp = [powerGAExp power];
proGAExp = [proGAExp propability];
scoreGAExp = [scoreGAExp fval];
% 提取最优解中的 phi 和 loc
optimalPhi = x(1:sat_num);
optimalidx = round(x(sat_num+1));
satSelection = selectionSatAll(optimalidx,:);
optimalXloc = sattochoose_x(satSelection);
optimalYloc = sattochoose_y(satSelection);
optimalSatTrack = sattrack_theta(satSelection);

%画出最优选星和相位情况下的功率分布
[aaa,corr] = power_calculation(sat_num,optimalXloc,optimalYloc,optimalPhi,x_ground,y_ground,H,node_x,node_y,D,frequency);
%powerdistribution(corr,x_ground,y_ground,node_x,node_y)

% 打印最终结果
disp(['最优解的相位值(pi)：', num2str(optimalPhi)])%pi
sat_locations = [];
for i =1:sat_num
sat_locations = [sat_locations optimalXloc(i) optimalYloc(i)];
end
disp(['最优解的卫星位置：',num2str(sat_locations)])
disp(['优化得分：',num2str(-fval)]) 
%对比无优化情况，最优选星，0相位作为基准
%x_base = [zeros(1,sat_num) optimalidx];
%[power_base] = singlenode_power(x_base,sat_num,x_ground,y_ground,node_x,node_y,sat_x,sat_y,frequency,D)

SATV = [sattochoose_x' sattochoose_y']/1000;
SATS = [optimalXloc' optimalYloc']/1000;
SATR = [sat_loc_xi' sat_loc_yi']/1000;
end
%%
%%%%得出选择的四颗卫星和各自的相位后，开始根据卫星运动实时更新卫星权重
sat_xmoving = v_sat*calculation_time*cosd(sattrack_theta(round(satSelection)));
sat_ymoving = v_sat*calculation_time*sind(sattrack_theta(round(satSelection)));
sat_xstart = optimalXloc;
sat_ystart = optimalYloc;%后续更新的初始卫星位置
phi = optimalPhi;%后续更新的初始相位

%算法更新参数初始化
time_update = 0.001;%更新间隔
visibleTime = 3;
times = 100;%更新次数
%times = visibleTime/time_update;

samplerate = 0.5;
test_update = samplerate*time_update;%监测真正的实时功率

%更高采样率的测试，用来检验无计算时的保持情况
test_t = times/samplerate;
sat_x_test = zeros(test_t,sat_num);
sat_y_test = zeros(test_t,sat_num);
power_test = zeros(test_t,1);
power_baseline = zeros(test_t,1);
propability_test = zeros(test_t,1);
propability_base = zeros(test_t,1);
score_test = zeros(test_t,1);
score_base = zeros(test_t,1);
for i = 1:test_t
sat_x_test(i,:) = sat_xstart+test_update*(i-1)*v_sat*cosd(sattrack_theta(round(satSelection)));
sat_y_test(i,:) = sat_ystart+test_update*(i-1)*v_sat*sind(sattrack_theta(round(satSelection)));
end
phi_test = zeros((1/samplerate)*times,sat_num);
phi_test(1,:) = wrapToPi(phi);
phi_base = zeros((1/samplerate)*times,sat_num);

%各项指标的初始化
score_record = [score_record fval];
power_record = [power_record power];
propability_record = [propability_record propability];

score_recordu0 = [score_recordu0 0];
propability_recordu0 = [propability_recordu0 0];
power_recordu0 = [power_recordu0 0];

score_recordGA = [score_recordGA fval];
power_recordGA = [power_recordGA power];
propability_recordGA = [propability_recordGA propability];
%各类变量的初始化
sat_x_record = zeros(times,sat_num);
sat_y_record = zeros(times,sat_num);%用来记录四颗卫星的轨迹
sat_x_record(1,:) = sat_xstart;
sat_y_record(1,:) = sat_ystart;
phi_record = zeros(times,sat_num);%用来记录四颗卫星的相位
%phi_record(1,:) = wrapToPi(phi);
phi_record(1,:) = phi;
phi_recordu0 = zeros(times,sat_num);%用来记录四颗卫星的相位
%phi_recordu0(1,:) = wrapToPi(phi);
phi_recordu0(1,:) = phi;
initialphi = ones(node_num,1)*phi_record(1,:);
%power_record = [power_record power];
%power_recordu0 = [power_recordu0 power];
R_now = [];
signalBF = [];
signalNoBF = [];
signalOriginal = [];

power_node_record = [];
poweru0_node_record = [];
powerGA_node_record = [];


for p =1:node_num
        for q= 1:sat_num
            R_now(p,q) = sqrt((x_ground(node_x(p))-sat_x_record(1,q)).^2 + (y_ground(node_y(p))-sat_y_record(1,q)).^2 + D^2);
        end
    end
    phi_temp = R_now*2*pi/lambda;
    signal_transmit = ones(node_num,sat_num);%卫星上的发射信号
    signal_transferred = signal_transmit.*exp(1i*phi_temp);%计算地面节点的接收信号(未叠加卫星上的发射相位，作为x)
    initGroundSignal = exp(1i*phi_record(1,:)) * signal_transferred.';%第一时刻的地面信号作为后续更新的理想信号




for n = 2:times

    %更新卫星位置
    sat_x_record(n,:) = sat_x_record(n-1,:) + v_sat*time_update*cosd(sattrack_theta(round(satSelection)));
    sat_y_record(n,:) = sat_y_record(n-1,:) + v_sat*time_update*sind(sattrack_theta(round(satSelection)));

    %计算路径长度
    for p =1:node_num
        for q= 1:sat_num
            R_now(p,q) = sqrt((x_ground(node_x(p))-sat_x_record(n,q)).^2 + (y_ground(node_y(p))-sat_y_record(n,q)).^2 + D^2);
        end
    end

    %计算路径相位
    phi_temp = R_now*2*pi/lambda;
    phi_SatTrack = 2*pi*optimalSatTrack/360;
    %计算地面节点的接收信号
    signal_transferred = signal_transmit.*exp(1i*phi_temp); 
    signal_receive = phi_record(n-1,:) * signal_transferred.';

    %根据当前卫星位置求解最优相位
    fun = @(x) objective(x,phi_temp,node_num,sat_num,power_base);
    x0 = phi_record(n-1,:);%以上一时刻的值作为起点
    options = optimoptions('fminunc','Display','off');%关闭计算状态显示
    tic
    [x,fval] = fminunc(fun,x0,options);%得出结果
    phi_record(n,:) = x;
    phi_recordu0(n,:) = phi_recordu0(n-1,:);%不做优化的对比情况
    signal_syn = exp(1i*phi_record(n,:)) * signal_transferred.';
    signal_synu0 = exp(1i*phi_recordu0(n,:)) * signal_transferred.';
    calculationmin = toc;%统计从计算开始到得出相位的时间
    calculationmin_record = [calculationmin_record calculationmin];
    %计算GA，作为算法 upper performance
    tic
    targetFunction = @(x1)adaptive_power(x1,sat_x_record(n,:),sat_y_record(n,:),sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D);
    nvars1 = sat_num;
    UB1 = pi*ones(1,sat_num);
    LB1 =-pi*ones(1,sat_num);
    options1 = optimoptions('ga','Display','off');
    [x1, fval] = ga(targetFunction, nvars1, [], [], [], [], LB1, UB1,[],options1);
    [fvalGA,powerGA,propabilityGA,powerGA_node] = adaptive_power(x1,sat_x_record(n,:),sat_y_record(n,:),sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D);
    calculationGA = toc;
    calculationGA_record = [calculationGA_record calculationGA];
    %记录GA结果
    power_recordGA = [power_recordGA powerGA];
    propability_recordGA = [propability_recordGA propabilityGA];
    score_recordGA = [score_recordGA fvalGA];
    
    %检验 如果算出来的值不如原始值 那就以原始值作为当前值
     if mean(signal_syn*signal_syn') < mean(signal_synu0*signal_synu0')
         phi_record(n,:) =  phi_recordu0(n,:);
     end
    

    %重新合成信号
    signal_syn = exp(1i*phi_record(n,:)) * signal_transferred.';
    signal_synu0 = exp(1i*phi_recordu0(n,:)) * signal_transferred.';
    
    %记录信号，与无优化状态做对比
    signalBF = [signalBF signal_syn.'];
    average_signalBF = mean(abs(signalBF),1);%求各个节点的平均幅度
    signalNoBF = [signalNoBF signal_synu0.'];
    average_signalNoBF = mean(abs(signalNoBF),1);
    signalOriginal = [signalOriginal sat_num*ones(1,node_num).'];

    %计算相关指标
    [score_temp, power_temp, propability_temp,power_node_temp]= adaptive_power(phi_record(n,:), sat_x_record(n,:), sat_y_record(n,:),sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D);
    [score_tempu0, power_tempu0, propability_tempu0,poweru0_node_temp]= adaptive_power(phi_recordu0(n,:), sat_x_record(n,:), sat_y_record(n,:),sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D);
    power_record = [power_record power_temp];
    power_recordu0 = [power_recordu0 power_tempu0];
    propability_record = [propability_record propability_temp];
    propability_recordu0 = [propability_recordu0 propability_tempu0];
    score_record = [score_record score_temp];
    score_recordu0 = [score_recordu0 score_tempu0];
    phi_test((n-1)/samplerate+1:n/samplerate,:) = x.*ones(1/samplerate,sat_num);%为更高采样率的画图做准备
    
    power_node_record = [power_node_record,power_node_temp];
    poweru0_node_record = [poweru0_node_record,poweru0_node_temp];
    powerGA_node_record = [powerGA_node_record,powerGA_node];
    %绘制卫星运动轨迹图
    %satelliteplot(sattochoose_x,sattochoose_y,sat_x_record,sat_y_record,n);

    %检验卫星是否仍然在可视范围内
    x_out = any(sat_x_record(n,:) > H_sat);
    y_out = any(sat_y_record(n,:) > H_sat);
    if (x_out) || (y_out)
        break%超出范围即结束当前选星
    end
end
%%
 power_node_record = 10*log10(power_node_record);
 poweru0_node_record = 10*log10(power_node_record);
 powerGA_node_record = 10*log10(power_node_record);
 power_record = 10*log10(power_record);
 power_recordu0 = 10*log10(power_recordu0);
 power_recordGA = 10*log10(power_recordGA);
%% 对节点概率做统计
data = power_node_record.';
% 绘制三维柱状图
figure;

b3 = bar3(data,0.7);
daspect([0.5 4 2]);
[x, y] = meshgrid(1:size(data, 2), 1:size(data, 1));
legend('Node.1','Node.2','Node.3','Node.4','Node.5','Node.6','Node.7','fontsize',30)
% 添加一个平行于 xy 平面的平面 z=4
% hold on;
% z_plane = ones(size(x)) * 4;
% surf(x, y, z_plane, 'FaceAlpha', 1, 'FaceColor', 'r', 'EdgeColor', 'none');
xlim([0 node_num+1])
ylim([1 times])
zlim([0 12])
 % 在 XZ 平面的 Y=39 处绘制一个平面
% hold on
% x_plane = [min(x(:)), min(x(:)); max(x(:)), max(x(:))]; % X=3
% y_plane = [39, 39; 39, 39]; % Y 的范围
% z_plane = [0, 16; 0, 16]; % Z 的范围
%surf(x_plane, y_plane, z_plane, 'FaceAlpha', 0.3); % 绘制平面，设置红色，透明度为0.3
set(gca, 'FontSize', 30);
xlabel('Node Label','fontsize', 50);
ylabel('Update Point','fontsize', 50);
zlabel('Power Gain(dB)','fontsize', 50);

%title('Power v.s Nodes and Update Points','fontsize', 50);

%%
figure
norMIMO = 10*log10(sat_num);
xMIMO = [0 length(data(39,:))+1 length(data(39,:))+1 0]; % 横坐标范围为 1 到 300
yMIMO = [0 0 norMIMO norMIMO];     % 纵坐标范围为 1 到 4
subplot(1,2,1)
b2 = bar(data(39,:),0.5,'FaceColor','flat');

set(gca, 'FontSize', 20);
%title('Node Power on One Update Point','fontsize', 30)
hold on
patch(xMIMO,yMIMO , 'r', 'FaceAlpha', 0.5); % 使用蓝色，设置半透明度为 0.5
xlabel('Node Label','fontsize', 50)
ylabel('Power Gain(dB)','fontsize', 50)
legend('','MIMO Baseline','fontsize', 30)
subplot(1,2,2)
xMIMO = [0 length(data(:,3))+1 length(data(:,3))+1 0];
bar(data(:,3),0.5,'cyan');
set(gca, 'FontSize', 30);
%title('Power of Node.3 v.s Update Points','fontsize', 30)
patch(xMIMO,yMIMO , 'r', 'FaceAlpha', 0.5);
xlabel('Update Point','fontsize', 50)
ylabel('Power Gain(dB)','fontsize', 50)
legend('','MIMO Baseline','fontsize', 30)

%%
% 统计小于 4 的个数
count4 = sum(data < 4, 2);

% 统计大于 4 但小于 8 的个数
count48 = sum(data >=4 & data < 8, 2);

% 统计大于 8 但小于 12 的个数
count812 = sum(data >= 8 & data < 12, 2);

% 统计大于 12 的个数
count12 = sum(data >= 12, 2);

count = [count4 count48 count812 count12];
data = count;
% 创建一个大小为 7x1 的逻辑向量，表示小于 4 的点
less_than_4 = any(data < 4, 2);

%% 计算积分功率
figure
power_integ = [];
poweru0_integ = [];
powerGA_integ = [];
    power_record = 10.^((power_record)/10);
    power_recordu0 = 10.^((power_recordu0)/10);
    power_recordGA = 10.^((power_recordGA)/10);
for i = 1:length(power_record)
    power_integ = [power_integ trapz(power_record(1:i))];
    poweru0_integ = [poweru0_integ trapz(power_recordu0(1:i))];
    powerGA_integ = [powerGA_integ trapz(power_recordGA(1:i))];
end

yyaxis left
plot(1:length(power_record),power_integ,'-gs',1:length(power_record),poweru0_integ,'-x',1:length(power_record),powerGA_integ,'-diamond','MarkerSize',7,'LineWidth',3)
yyaxis right
bar(power_record(1:i),0.7);

yyaxis left
ylabel('Receiving Energy','fontsize',50)
yyaxis right
ylabel('Normalized Power','fontsize',50)
hold on
x = [0 300 300 0]; % 横坐标范围为 1 到 300
y = [0 0 4 4];     % 纵坐标范围为 1 到 4
% 绘制半透明平面
patch(x, y, 'blue', 'FaceAlpha', 0.5); % 使用蓝色，设置半透明度为 0.5
set(gca, 'fontsize', 30)
grid on
xlim([0,300])
ylim([0,16])
legend('Proposed Algorithm','No Optimization','GA','fontsize',30)
xlabel('Update Points','fontsize',50)
%title('Average Receiving Energy','fontsize',40)
%%
    %画出计算时间对比曲线
    plot(1:length(calculationGA_record),calculationmin_record,1:length(calculationGA_record),calculationGA_record,'LineWidth',3)
    title('计算时间对比')
    legend('minlinear','GA')
%%
    %显示平均值
    disp('整个时段内')
    disp(['平均功率：',num2str(mean(power_record))])
    disp(['平均概率：',num2str(mean(propability_record))]) 
    disp(['平均分数：',num2str(mean(-score_record))])
    disp(['平均功率GA：',num2str(mean(power_recordGA))])
    disp(['平均概率GA：',num2str(mean(propability_recordGA))]) 
    disp(['平均分数GA：',num2str(mean(-score_recordGA))])
         
% 计算更高采样率下的各项参数
for i = 1:test_t
[score_test(i),power_test(i),propability_test(i)] = adaptive_power(phi_test(i,:),sat_x_test(i,:),sat_y_test(i,:),sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D);
[score_base(i),power_baseline(i),propability_base(i)] = adaptive_power(phi_base(i,:),sat_x_test(i,:),sat_y_test(i,:),sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D);
end

%在最后一次选星中画出实时参数对比
if a == selectMax
figure
plot(1:n/samplerate,power_test(1:n/samplerate),1:n/samplerate,power_baseline(1:n/samplerate), 'LineWidth', 2)%在超出范围后截断
title('实时功率')
legend('优化实时功率','无优化实时功率')
figure
plot(1:n/samplerate,propability_test(1:n/samplerate),1:n/samplerate,propability_base(1:n/samplerate), 'LineWidth', 2)%在超出范围后截断
title('实时概率')
legend('优化实时概率','无优化实时概率')
figure
plot(1:n/samplerate,-score_test(1:n/samplerate),1:n/samplerate,-score_base(1:n/samplerate), 'LineWidth', 2)%在超出范围后截断
title('实时分数')
legend('优化实时分数','无优化实时分数')

%各个节点的接收信号强度对比，比较五优化和优化后的结果
figure
for i = 1:node_num
subplot(1,node_num,i)
plot(1:size(signalBF(i,:),2),abs(signalBF(i,:)),'b','LineWidth',1,'LineStyle','--');
hold on 
plot(1:size(signalOriginal(i,:),2),abs(signalOriginal(i,:)),'r','LineWidth',1,'LineStyle',':');
hold on 
plot(1:size(signalNoBF(i,:),2),abs(signalNoBF(i,:)),'LineWidth',1,'LineStyle','-.');
legend('signalBF','signalOriginal','signalNoBF')
end
title('各个节点接收信号')

%节点接收信号的平均值对比
figure
plot(1:length(average_signalBF),average_signalBF,'LineWidth',1)
hold on
plot(1:length(average_signalNoBF),average_signalNoBF,'LineWidth',1)
legend('avergaeBF','averageNoBF')
title('节点平均对比')

end

%记录每次选定卫星后持续了多久才进行的下一次选星
outtime_record = [outtime_record n];

% 选星循环
%end 

%%
% 画出各项参数的分段曲线图，描述整个过程内各项指标的变化情况，baseline为无优化下的各项指标

plot_segmented(outtime_record,power_record,power_recordu0);%绘制功率的分段曲线
title('Power Record Segments');
xlabel('Data Point');
ylabel('Power');

plot_segmented(outtime_record,propability_record,propability_recordu0);%绘制功率的分段曲线
title(' Propability Record Segments');
xlabel('Data Point');
ylabel('Power');

plot_segmented(outtime_record,-score_record,-score_recordu0);%绘制功率的分段曲线
title('Score Record Segments');
xlabel('Data Point');
ylabel('Power');

plot_segmented(outtime_record,power_record,power_recordGA);%绘制功率的分段曲线
title('Power minlinear vs GA');
xlabel('Data Point');
ylabel('Power');

plot_segmented(outtime_record,propability_record,propability_recordGA);%绘制功率的分段曲线
title(' Propability minlinear vs GA');
xlabel('Data Point');
ylabel('Power');

plot_segmented(outtime_record,-score_record,-score_recordGA);%绘制功率的分段曲线
title('Score minlinear vs GA');
xlabel('Data Point');
ylabel('Power');
















%%%%%%%%%%%%%%%%%%%%%%%函数部分
%%%目标函数，先对路径和卫星初始相位求和并归一化到正负Pi，再求均值误差
function result = objective(x, phi_temp,node_num,sat_num,power_base)
    result = 0;
    pow = 0;
    pro = 0;
    pro_sigmoid = 0;
    for i = 1:node_num
    %归一化到【-pi pi】
    %phi_normalize = wrapToPi(x + phi_temp(i,:));
    %difNode = abs(phi_normalize-mean(phi_normalize));%每个节点四条路径间的相位差异
    %difNode = sum(difNode);%求出每个节点的总差异
    %result = result + difNode;%叠加所有节点

    %映射到exp
    difNode = exp(1i*(x + phi_temp(i,:)));
    difNode = abs(sum(difNode));%求出每个节点的总差异
    if difNode > sqrt(power_base)
        pro = pro + 1;
    end
    pro_sigmoid = pro_sigmoid + 1/(1+exp(-difNode+sat_num/2));%将节点信号强度调整到中间值为0的分布，在映射到sigmoid函数
    pow = pow - difNode;%叠加所有节点(求最小故取负)
    end
    result = 0.5 * (pow/(sat_num*node_num)) - 0.5 * pro_sigmoid/node_num;
    %result = pow;
end


%%%单节点功率计算
function [power] = singlenode_power(xx,sat_num,x_ground,y_ground,node_x,node_y,x_sat,y_sat,frequency,D)
phi = xx(1:sat_num);
satidx = xx(sat_num+1:sat_num*2);  
sat_loc_x = x_sat(round(satidx));
sat_loc_y = y_sat(round(satidx));

A = ones(1,sat_num);
%satellite coordinates 
c = 3e8;
lambda = c/frequency;%波长
R = zeros(1,sat_num);
delta_phi = 0;
signal_receive = 0;
%%
   for m = 1:length(node_x)
        for k = 1:sat_num
            R(k) = sqrt((y_ground(node_y(m))-sat_loc_y(k)).^2+(x_ground(node_x(m))-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            signal_receive = signal_receive + A(k).*exp(1i*(delta_phi(k)+phi(k)));
        end
        cor(m) = signal_receive*signal_receive';
        signal_receive = 0;
   end
   power = mean(cor);
end


%%%各项指标统计
function [score,power,propability] = score(xx,sat_num,x_ground,y_ground,node_x,node_y,x_sat,y_sat,power_base,frequency,D,selectionSatAll)

phi = xx(1:sat_num);
satidx = round(xx(sat_num+1));  
satSelection = selectionSatAll(satidx,:);
sat_loc_x = x_sat(satSelection);
sat_loc_y = y_sat(satSelection);

A = ones(1,sat_num);
c = 3e8;
lambda = c/frequency;%波长
R = zeros(1,sat_num);
delta_phi = 0;
signal_receive = 0;
propability = 0;
propability_sigmoid = 0;

   for m = 1:length(node_x)
        for k = 1:sat_num
            R(k) = sqrt((y_ground(node_y(m))-sat_loc_y(k)).^2+(x_ground(node_x(m))-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            signal_receive = signal_receive + A(k).*exp(1i*(delta_phi(k)+phi(k)));
        end
        cor(m) = signal_receive*signal_receive';

        signal_receive = 0;
        if cor(m)>power_base
            propability = propability + 1;
        end
        propability_sigmoid = propability_sigmoid + 1/(1+exp(-cor(m)+sat_num^2/2));
   end
   propability = (propability)/length(node_x);
   propability_sigmoid = (propability_sigmoid)/length(node_x);
   power = mean(cor);
   score = -0.5*power/sat_num.^2-0.5*propability_sigmoid;
end


%%% 分段指标绘制函数
function plot_segmented(outtime_record,record,base)
% Initialize figure
figure;
startIndex = 1;  % Start index for each segment
for i = 1:length(outtime_record)
    endIndex = startIndex + outtime_record(i) - 1;  % End index for each segment
    % Extract segment from power_record
    segment = record(startIndex:endIndex);
    % Plot the segment
    plot(startIndex:endIndex, segment, 'LineWidth', 2,'Color','b');
    hold on;
    baseline = base(startIndex:endIndex);%基准线，无优化时的功率
    plot(startIndex:endIndex, baseline, 'LineWidth', 2,'Color','r');
    y = 0:0.1:max(record);
    plot(endIndex*ones(size(y)),y,'LineWidth', 2,'Color','g');
    hold on;  % Hold the plot for the next segment
    startIndex = endIndex + 1;  % Update the start index for the next segment
end
hold off;  % Release the hold after all segments are plotted
end


%%% 卫星动态情况绘制函数
function satelliteplot(sattochoose_x,sattochoose_y,sat_x_record,sat_y_record,n)
    scatter(sattochoose_x,sattochoose_y,'x')
    hold on
    scatter(sat_x_record(n,:), sat_y_record(n,:),50,'b');
    legend('所有可选卫星','通过算法选取')
    pause(0.5)
    drawnow
    clf
end


%%% 与score函数的差别在于输入x只有相位，长度为sat_num
function [score,power,propability,power_node] = adaptive_power(x,sat_loc_x,sat_loc_y,sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D)

phi = x;
A = ones(1,sat_num);
c = 3e8;
lambda = c/frequency;%波长
R = zeros(1,sat_num);
delta_phi = 0;
signal_receive = 0;
propability = 0;
propability_sigmoid = 0;
power_node = zeros(length(node_x),1);
   for m = 1:length(node_x)
        for k = 1:sat_num
            R(k) = sqrt((y_ground(node_y(m))-sat_loc_y(k)).^2+(x_ground(node_x(m))-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            signal_receive = signal_receive + A(k).*exp(1i*(delta_phi(k)+phi(k)));
        end
        cor(m) = signal_receive*signal_receive';
        power_node(m) = cor(m);
        signal_receive = 0;
        cor_normalized = cor(m)-sat_num^2/2;
        propability_sigmoid = propability_sigmoid + 1/(1+exp(-cor_normalized));
         if cor(m)>power_base
            propability = propability + 1;
         end
   end
   propability = (propability)/length(node_x);
   propability_sigmoid = (propability_sigmoid)/length(node_x);
   power = mean(cor);
   score = -0.5*power/(sat_num^2)-0.5*propability_sigmoid;
end