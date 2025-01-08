clc
clear
close all
% 优化的变量：选星（从给定的sat_available = 10个卫星中选取sat_num = 4个卫星）； 每颗卫星的相位phi
% 解的构成：[phi*sat_num index*sat_num(选取的卫星标号)]
% 目标函数：地面端单个节点的功率，考虑传输过程中的衰减
% 优化算法：ga/pso
% 选星区域：直径H_sat = 10km的范围内 
% 轨道高度 D = 550km

% 计算选星和相位
% 根据轨道计算卫星位置和相位
% 补偿该相位变化
% basic parameters

c = 3e8;
frequency = 915e6;
H_sat = 10; %卫星分布范围
H = 20;%地面范围
w = frequency*2*pi;
lambda = c/frequency;%波长
%lambda = lambda/1000;%换算成km
sat_num = 4;
v_sat = 1;%卫星运动速度为7.8km/s
D = 50;%orbit 
node_num = 1;

step = H/201;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;

node_x = randi(length(x_ground),1,node_num);

node_y = randi(length(y_ground),1,node_num);


step_sat = H_sat/201;
sat_available = 10;

sat_x = -H_sat/2:step_sat:H_sat/2;
sat_y = -H_sat/2:step_sat:H_sat/2;

nvars = 2*sat_num;%优化目标为相位和四颗卫星的标号
power_record = [];%用来记录每次更新的地面功率
power_recordu0 = [];%用来记录权重u=0时的baseline
outtime_record = [];%用来记录每次选星后的持续时间

figure%卫星图
timeupdate = [];
uu = [];
for a = 1:1
% 产生可选卫星的坐标
sattochoose_x = sat_x(randperm(length(sat_x), sat_available));%可选卫星的坐标
sattochoose_y = sat_y(randperm(length(sat_y), sat_available));
sattrack_theta =  0 + rand(1,sat_available) * (90 - 0);%0到90度内的随机数
estimate_calculation_time = 0;%初始化计算时间/(以估计的运算时间及逆行计算)
% 使用 ga 函数进行单目标优化
%随即选星，0相位作为基准
x_base = [zeros(1,sat_num) randi(sat_available,1,sat_num)];
[power_base] = singlenode_power(x_base,sat_num,x_ground,y_ground,node_x,node_y,sattochoose_x,sattochoose_y,frequency,D)
phii = x_base(1:sat_num);
satidx = x_base(sat_num+1:sat_num*2);  
sat_loc_xi = sattochoose_x(round(satidx));
sat_loc_yi = sattochoose_y(round(satidx));
[aaa,corr] = power_calculation(sat_num,sat_loc_xi,sat_loc_yi,phii,x_ground,y_ground,H,node_x,node_y,D,frequency);
powerdistribution(corr,x_ground,y_ground,node_x,node_y)
%%
targetFunction = @(x)score(x,sat_num,x_ground,y_ground,node_x,node_y,sattochoose_x,sattochoose_y,power_base,frequency,D);
%为了避免重复，将选星分为四个区间，每个变量从一个区间中选取
UB =[pi*ones(1,sat_num)   round(sat_available/4) round(sat_available/2)      round(sat_available*3/4)    sat_available];
LB =[-pi*ones(1,sat_num)  1                      round(sat_available/4)+0.01 round(sat_available/2)+0.01 round(sat_available*3/4)+0.01];

%options = optimoptions('ga','PopulationSize',50,'MaxGenerations',10);%default: ps 200; mg 200*number of variable
options = optimoptions('ga','PopulationSize',100);
tic;
[x, fval,exitflag,output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,[],options);
calculation_time = toc
[fval,power,propability] = score(x,sat_num,x_ground,y_ground,node_x,node_y,sattochoose_x,sattochoose_y,power_base,frequency,D);

% 提取最优解中的 phi 和 loc
optimalPhi = x(1:sat_num);
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
disp(['优化得分：',num2str(-fval)]) 
%对比无优化情况，最优选星，0相位作为基准
%x_base = [zeros(1,sat_num) optimalidx];
%[power_base] = singlenode_power(x_base,sat_num,x_ground,y_ground,node_x,node_y,sat_x,sat_y,frequency,D)
disp(['无优化的功率值（选星相同，但不改变相位）：',num2str(power_base)])


%%%%得出选择的四颗卫星和各自的相位后，开始根据卫星运动实时更新卫星权重
sat_xmoving = v_sat*calculation_time*cosd(sattrack_theta(round(optimalidx)));
sat_ymoving = v_sat*calculation_time*sind(sattrack_theta(round(optimalidx)));
sat_xstart = optimalXloc;%+sat_xmoving;
sat_ystart = optimalYloc;%+sat_ymoving;%后续更新的初始卫星位置
phi = optimalPhi;%后续更新的初始相位
%phi = mod(phi,2*pi);
%ep = zeros(10,10);
% for aa = 1:100
%     for bb  = 1:100
%         power_record = [];
%     time_update = 0.001*aa;%更新间隔
%     times = 100;%更新次数
%     u = 0.001*bb;%更新权重    

time_update = 0.01;%更新间隔
times = 1000;%更新次数

sat_x_record = zeros(times,sat_num);
sat_y_record = zeros(times,sat_num);%用来记录四颗卫星的轨迹
sat_x_record(1,:) = sat_xstart;
sat_y_record(1,:) = sat_ystart;
phi_record = zeros(times,sat_num);%用来记录四颗卫星的相位
%phi_record(1,:) = wrapToPi(phi);
phi_record(1,:) = exp(1i*phi);
phi_recordu0 = zeros(times,sat_num);%用来记录四颗卫星的相位
%phi_recordu0(1,:) = wrapToPi(phi);
phi_recordu0(1,:) = exp(1i*phi);
initialphi = ones(node_num,1)*phi_record(1,:);
%power_record = [power_record power];
%power_recordu0 = [power_recordu0 power];
error_record = [];
R_now = [];
dif = [];%记录信号差异
signalBF = [];
signalNoBF = [];
signalOriginal = [];
%MI_record = [];
%MO_record = [];
%LMS自适应更新phi，接收信号为地面节点的信号，参考信号为忽略路程相移的信号
%n=1时刻为GA进行时的卫星位置和计算出来的相位
%n=2时刻
sum_error = 0;%error累计初始化

for p =1:node_num
        for q= 1:sat_num
            R_now(p,q) = sqrt((x_ground(node_x(p))-sat_x_record(1,q)).^2 + (y_ground(node_y(p))-sat_y_record(1,q)).^2 + D^2);
        end
    end
    phi_temp = exp(1i*R_now*2*pi/lambda);
    signal_transmit = ones(node_num,sat_num);%卫星上的发射信号
    signal_transferred = signal_transmit.*phi_temp;%计算地面节点的接收信号(未叠加卫星上的发射相位，作为x)
    initGroundSignal = phi_record(1,:) * signal_transferred.';%第一时刻的地面信号作为后续更新的理想信号

for n = 2:times
    
    sat_x_record(n,:) = sat_x_record(n-1,:) + v_sat*time_update*cosd(sattrack_theta(round(optimalidx)));
    sat_y_record(n,:) = sat_y_record(n-1,:) + v_sat*time_update*sind(sattrack_theta(round(optimalidx)));
   
    
    for p =1:node_num
        for q= 1:sat_num
            R_now(p,q) = sqrt((x_ground(node_x(p))-sat_x_record(n,q)).^2 + (y_ground(node_y(p))-sat_y_record(n,q)).^2 + D^2);
        end
    end
    phi_temp = exp(1i*R_now*2*pi/lambda);

    % H = exp(1j*(phi_temp));
    % [U,S,V] = svd(H);
    % signal_MI = ones(1,sat_num)*exp(1i*w*n*time_update);
    % MI_record = [MI_record signal_MI];
    % signal_MI = V*signal_MI';
    % signal_MO = U^-1 * H* signal_MI;
    % signal_MO = signal_MO./(diag(S));
    % MO_record = [MO_record signal_MO];
    


    signal_transferred = signal_transmit.*phi_temp %计算地面节点的接收信号(未叠加卫星上的发射相位，作为x)
    signal_receive = phi_record(n-1,:) * signal_transferred.';
  
    %梯度下降
    % for k=1:node_num
    %     grad(k,:) = signal_transferred(k,:).*(signal_receive(k)-initGroundSignal*exp(1i*w*(n-1)*time_update));
    % end
    % gradallnode = sum(grad,1);
    % u = 0.01;
    % u0 = 0;
    % phi_record(n,:) = phi_record(n-1,:) - 1/2*u*gradallnode;
    % phi_record(n,:) = exp(1i*angle(phi_record(n,:)));
    % phi_recordu0(n,:) = phi_recordu0(n-1,:);
    


    %自适应LMS 
    d_r = initGroundSignal;%期望响应
    error = d_r - sum(abs(signal_receive))/node_num;%对信号取abs是因为信号传输路径会产生相位差
    error_record = [error_record error]; 
    % % %更新权重u(LMS)
    % average_error = mean(error_record);
    % f = exp(average_error)/(exp(average_error)+1);
    % sum_error = sum_error + (abs(error));
    % node_signal = signal_transferred;%各个节点的合成信号
    % estimated_power = node_signal*node_signal';%各节点平均功率
    % u = f/(sum_error^2+estimated_power);

    miu = max(diag(signal_receive*signal_receive'));
    u = 0.01;
    phi_record(n,:) = phi_record(n-1,:) + u*error.'*signal_transferred;
    phi_record(n,:) = exp(1i*angle(phi_record(n,:)));
    u0 = 0;
    phi_recordu0(n,:) = phi_recordu0(n-1,:) + u0*error.'*signal_transferred;
    
    %归一化到【-pi pi】, 求功率
    %phi_record(n,:) = wrapToPi(phi_record(n,:));
    %phi_recordu0(n,:) =  wrapToPi(phi_recordu0(n,:));
    signalBF = [signalBF phi_record(n,:) * signal_transferred.'];
    signalNoBF = [signalNoBF phi_recordu0(n,:) * signal_transferred.'];
    signalOriginal = [signalOriginal initGroundSignal];
    dif = [dif phi_record(n,:) * signal_transferred.' - initGroundSignal];
    [score_temp, power_temp, propability_temp]= adaptive_power(phi_record(n,:), sat_x_record(n,:), sat_y_record(n,:),sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D);
    [score_tempu0, power_tempu0, propability_tempu0]= adaptive_power(phi_recordu0(n,:), sat_x_record(n,:), sat_y_record(n,:),sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D);
    power_record = [power_record power_temp];
    power_recordu0 = [power_recordu0 power_tempu0];
    x_out = any(sat_x_record(n,:) > H_sat);
    y_out = any(sat_y_record(n,:) > H_sat);
    %satelliteplot(sattochoose_x,sattochoose_y,sat_x_record,sat_y_record,n);%绘制卫星运动轨迹图
    %ep(aa,bb) = sum(abs(power_record - power_base))/length(power_record);
    

    if (x_out) || (y_out)
        break
    end
end
%% MIMO figure
% figure
% for i = 1:size(MO_record,1)
% plot(1:size(MO_record,2),abs(MO_record(i,:)))
% hold on
% end
% hold on
% plot(1:size(MI_record,1),abs(MI_record))
% legend('signal1','2','3','4','MI')
%% error figure
figure
subplot(1,2,1)
plot(1:size(phi_record,1),angle(phi_record(:,1)));
title('phi')
subplot(1,2,2)
plot(1:length(dif),dif);
title('dif')

figure
subplot(1,2,1)
plot(1:size(signalBF,2),abs(signalBF));
title('signalBF')
hold on 
plot(1:size(signalOriginal,2),abs(signalOriginal));
subplot(1,2,2)
plot(1:size(signalNoBF,2),abs(signalNoBF));
hold on 
plot(1:size(signalOriginal,2),abs(signalOriginal));
title('signalNoBF')
%     end
% end
outtime_record = [outtime_record n];%即每次选定卫星后持续了多久再进行下一次选星
% value = min(min(ep));
% [row col] = find(ep==value);
% timeupdate = [timeupdate row*0.001];
% uu = [uu col*0.001];
end
%powerplot_segmented(outtime_record,power_record,power_recordu0,-power_base);%绘制功率的分段曲线






















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

%satellite coordinates 


RNM = zeros(1,sat_num);
R = zeros(1,sat_num);
phi_moving = zeros(1,sat_num);
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


function [score,power,propability] = score(xx,sat_num,x_ground,y_ground,node_x,node_y,x_sat,y_sat,power_base,frequency,D)
phi = xx(1:sat_num);
satidx = xx(sat_num+1:sat_num*2);  
sat_loc_x = x_sat(round(satidx));
sat_loc_y = y_sat(round(satidx));
A = ones(1,sat_num);
%satellite coordinates 
c = 3e8;
lambda = c/frequency;%波长

%satellite coordinates 

RNM = zeros(1,sat_num);
R = zeros(1,sat_num);
phi_moving = zeros(1,sat_num);
delta_phi = 0;
signal_receive = 0;
propability = 0;
%%
   for m = 1:length(node_x)
        for k = 1:sat_num
            R(k) = sqrt((y_ground(node_y(m))-sat_loc_y(k)).^2+(x_ground(node_x(m))-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            signal_receive = signal_receive + A(k).*exp(1i*(delta_phi(k)+phi(k)));
        end
        cor(m) = signal_receive*signal_receive';
        signal_receive = 0;
        if cor(m)>power_base
            propability = propability+1;
        end
   end
   propability = (propability)/length(node_x);
   power = mean(cor);
   score = -0.5*power/sat_num.^2-0.5*propability;
end


function powerplot_segmented(outtime_record,power_record,power_recordu0,power_base)

% Initialize figure
figure;

startIndex = 1;  % Start index for each segment

for i = 1:length(outtime_record)
    endIndex = startIndex + outtime_record(i) - 1;  % End index for each segment
    
    % Extract segment from power_record
    segment_power = power_record(startIndex:endIndex);
    segment_poweru0 = power_recordu0(startIndex:endIndex);
    % Plot the segment
    plot(startIndex:endIndex, segment_power, 'LineWidth', 2,'Color','b');
    hold on;
    plot(startIndex:endIndex, segment_poweru0, 'LineWidth', 2,'Color','r');
    hold on;
    baseline = power_base*ones(1,length(startIndex:endIndex));%基准线，无优化时的功率
    plot(startIndex:endIndex, -baseline, 'LineWidth', 2);
    hold on;  % Hold the plot for the next segment
    
    startIndex = endIndex + 1;  % Update the start index for the next segment
end

hold off;  % Release the hold after all segments are plotted

title('Power Record Segments');
xlabel('Data Point');
ylabel('Power');

end

function satelliteplot(sattochoose_x,sattochoose_y,sat_x_record,sat_y_record,n)
    scatter(sattochoose_x,sattochoose_y,'x')
    hold on
    scatter(sat_x_record(n,:), sat_y_record(n,:),50,'b');
    legend('所有可选卫星','通过算法选取')
    pause(0.5)
    drawnow
    clf
end

function [score,power,propability] = adaptive_power(x,sat_loc_x,sat_loc_y,sat_num,x_ground,y_ground,node_x,node_y,power_base,frequency,D)
phi = wrapToPi(angle(x));
%phi = x;
A = abs(x);
%A = ones(1,sat_num);

%satellite coordinates 
c = 3e8;
lambda = c/frequency;%波长
%satellite coordinates 
RNM = zeros(1,sat_num);
R = zeros(1,sat_num);
phi_moving = zeros(1,sat_num);
delta_phi = 0;
signal_receive = 0;
propability = 0;
%%
   for m = 1:length(node_x)
        for k = 1:sat_num
            R(k) = sqrt((y_ground(node_y(m))-sat_loc_y(k)).^2+(x_ground(node_x(m))-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            signal_receive = signal_receive + A(k).*exp(1i*(delta_phi(k)+phi(k)));
        end
        bb = signal_receive;
        cor(m) = signal_receive*signal_receive';
        signal_receive = 0;
         if cor(m)>power_base
            propability = propability+1;
         end
     
   end
   
   propability = (propability)/length(node_x);
   power = mean(cor);
   score = -0.5*power/sat_num.^2-0.5*propability;
end