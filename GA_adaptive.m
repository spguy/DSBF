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
frequency = 3.5e9;
w = frequency*2*pi;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km

sat_num = 4;

H = 1;%地面范围
step = H/2001;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;
node_num = 10;
node_x = randi(length(x_ground),1,node_num);
%node_x = [488	343	520	1837	1530	547	1323	150	1358	346];
node_y = randi(length(y_ground),1,node_num);
%node_y = [776	1609	758	1651	280	303	1921	1757	1732	625]; 
H_sat = 10; %卫星分布范围
step_sat = H_sat/201;
sat_available = 10;
v_sat = 7.8;%卫星运动速度为7.8km/s
sat_x = -H_sat/2:step_sat:H_sat/2;
sat_y = -H_sat/2:step_sat:H_sat/2;
D = 550;%orbit 
nvars = 2*sat_num;%优化目标为相位和四颗卫星的标号
power_record = [];%用来记录每次更新的地面功率
score_record = [];
propability_record = [];
outtime_record = [];%用来记录每次选星后的持续时间

figure%卫星图
timeupdate = [];
uu = [];
for a = 1:5
% 产生可选卫星的坐标
sattochoose_x = sat_x(randperm(length(sat_x), sat_available));%可选卫星的坐标
sattochoose_y = sat_y(randperm(length(sat_y), sat_available));
%sattochoose_x = [-3.95522388059702	0.422885572139303	3.20895522388060	0.820895522388059	0.273631840796019	0.124378109452736	-4.60199004975124	-0.572139303482587	4.15422885572139	4.45273631840796];
%sattochoose_y = [3.05970149253731	-1.31840796019900	0.621890547263681	-4.20398009950249	-3.60696517412935	0.721393034825870	2.06467661691542	2.86069651741294	-4.80099502487562	3.15920398009950]
sattrack_theta =  0 + rand(1,sat_available) * (90 - 0);%0到90度内的随机数
%sattrack_theta = [24.3095282693811	18.4963733898125	4.17218957452401	1.99620145982118	65.3223917838014	67.4720062288278	30.3667178141722	34.5312650582861	74.6880294400745	300000000];

estimate_calculation_time = 0;%初始化计算时间/(以估计的运算时间及逆行计算)
calculate_limit = 0.01;%最大计算时间
% 使用 ga 函数进行单目标优化
x_base = [zeros(1,sat_num) randi(sat_available,1,sat_num)];
[power_base] = singlenode_power(x_base,sat_num,x_ground,y_ground,node_x,node_y,sat_x,sat_y);
targetFunction = @(x)score(x,sat_num,x_ground,y_ground,node_x,node_y,sattochoose_x,sattochoose_y,power_base);
%为了避免重复，将选星分为四个区间，每个变量从一个区间中选取
UB =[pi*ones(1,sat_num)   round(sat_available/4) round(sat_available/2)      round(sat_available*3/4)    sat_available];
LB =[-pi*ones(1,sat_num)  1                      round(sat_available/4)+0.01 round(sat_available/2)+0.01 round(sat_available*3/4)+0.01];

%options = optimoptions('ga','PopulationSize',50,'MaxGenerations',10);%default: ps 200; mg 200*number of variable
options = optimoptions('ga','PopulationSize',100,'MaxTime',calculate_limit);
tic;
[x, fval,exitflag,output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,[],options);
calculation_time = toc
[fval,power,propability] = score(x,sat_num,x_ground,y_ground,node_x,node_y,sattochoose_x,sattochoose_y,power_base);


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
disp(['节点处功率值(NM)：',num2str(power)]) 
disp(['超过门限概率：',num2str(propability)])


disp(['无优化的功率值（选星相同，但不改变相位）：',num2str(power_base)])
%%
% [power_moving,sat_xmoving,sat_ymoving,phi_moving] = singlenode_power(x,sat_num,v_sat,calculation_time,sattrack_theta,node_x,node_y,x_sat,y_sat);
% disp(['节点处功率值(M)：',num2str(-power_moving)])
% disp(['卫星x方向位移',num2str(sat_xmoving)])
% disp(['卫星y方向位移',num2str(sat_ymoving)])
% disp(['卫星位移',num2str(calculation_time*v_sat)])
% disp(['相位变化',num2str(phi_moving)])

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

time_update = calculate_limit;%更新间隔
times = 500;%更新次数
samplerate = 0.5;
test_update = samplerate*time_update;%监测真正的实时功率

sat_x_record = zeros(times,sat_num);
sat_y_record = zeros(times,sat_num);%用来记录四颗卫星的轨迹
sat_x_record(1,:) = sat_xstart;
sat_y_record(1,:) = sat_ystart;
phi_record = zeros(times,sat_num);%用来记录四颗卫星的相位
phi_record(1,:) = wrapToPi(phi);

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
sat_x_test(i,:) = sat_xstart+test_update*(i-1)*v_sat*cosd(sattrack_theta(round(optimalidx)));
sat_y_test(i,:) = sat_ystart+test_update*(i-1)*v_sat*sind(sattrack_theta(round(optimalidx)));
end
phi_test = zeros((1/samplerate)*times,sat_num);%用来记录四颗卫星的相位
phi_test(1,:) = wrapToPi(phi);
phi_base = zeros((1/samplerate)*times,sat_num);

initialphi = ones(node_num,1)*phi_record(1,:);
score_record = [score_record fval];
power_record = [power_record power];
propability_record = [propability_record propability];
R_now = [];
%LMS自适应更新phi，接收信号为地面节点的信号，参考信号为忽略路程相移的信号
%n=1时刻为GA进行时的卫星位置和计算出来的相位
%n=2时刻
sum_error = 0;%error累计初始化

for n = 2:times
    signal_transmit = ones(node_num,sat_num).*exp(1j*w*(n-1)*time_update);%卫星上的发射信号
    sat_x_record(n,:) = sat_x_record(n-1,:) + v_sat*time_update*cosd(sattrack_theta(round(optimalidx)));
    sat_y_record(n,:) = sat_y_record(n-1,:) + v_sat*time_update*sind(sattrack_theta(round(optimalidx)));
    targetFunction = @(x1)adaptive_power(x1,sat_x_record(n,:),sat_y_record(n,:),sat_num,x_ground,y_ground,node_x,node_y,power_base);
    nvars1 = sat_num;
    UB1 = pi*ones(1,sat_num);
    LB1 =-pi*ones(1,sat_num);
    options1 = optimoptions('ga','MaxTime',calculate_limit);
    [x1, fval] = ga(targetFunction, nvars1, [], [], [], [], LB1, UB1,[],options1);
    [fval,power,propability] = adaptive_power(x1,sat_x_record(n,:),sat_y_record(n,:),sat_num,x_ground,y_ground,node_x,node_y,power_base);
    score_record = [score_record fval];
    power_record = [power_record power];%代表每次计算后最优功率的变化
    propability_record = [propability_record propability];
    phi_test((n-1)/samplerate+1:n/samplerate,:) = x1.*ones(1/samplerate,sat_num);
    x_out = any(sat_x_record(n,:) > H_sat);
    y_out = any(sat_y_record(n,:) > H_sat);
    %satelliteplot(sattochoose_x,sattochoose_y,sat_x_record,sat_y_record,n);%绘制卫星运动轨迹图
    %ep(aa,bb) = sum(abs(power_record - power_base))/length(power_record);
    if (x_out) || (y_out)
        break
    end
end
for i = 1:test_t
[score_test(i),power_test(i),propability_test(i)] = adaptive_power(phi_test(i,:),sat_x_test(i,:),sat_y_test(i,:),sat_num,x_ground,y_ground,node_x,node_y,power_base);
[score_base(i),power_baseline(i),propability_base(i)] = adaptive_power(phi_base(i,:),sat_x_test(i,:),sat_y_test(i,:),sat_num,x_ground,y_ground,node_x,node_y,power_base);
end
%%
if a == 5
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
end
% end
% end
outtime_record = [outtime_record n];%即每次选定卫星后持续了多久再进行下一次选星
% value = min(min(ep));
% [row col] = find(ep==value);
% timeupdate = [timeupdate row*0.001];
% uu = [uu col*0.001];
end
%%
plot_segmented(outtime_record,power_record,power_baseline);%绘制功率的分段曲线
title('Power Record Segments');
xlabel('Data Point');
ylabel('Power');

plot_segmented(outtime_record,propability_record,propability_base);%绘制功率的分段曲线
title(' Propability Record Segments');
xlabel('Data Point');
ylabel('Power');

plot_segmented(outtime_record,-score_record,-score_base);%绘制功率的分段曲线
title('Score Record Segments');
xlabel('Data Point');
ylabel('Power');






















%%%单节点功率计算
function [power,sat_xmoving,sat_ymoving,phi_moving] = singlenode_power(xx,sat_num,x_ground,y_ground,node_x,node_y,x_sat,y_sat)
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
D = 550;
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


function [score,power,propability] = score(xx,sat_num,x_ground,y_ground,node_x,node_y,x_sat,y_sat,power_base)
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
D = 550;
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


%%动态功率计算
function [score,power,propability] = adaptive_power(x,sat_loc_x,sat_loc_y,sat_num,x_ground,y_ground,node_x,node_y,power_base)
phi = x;
A = ones(1,sat_num);
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
D = 550;
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


%更新运动过程的power
function [power] = powerrecording(sat_num,x_ground,y_ground,node_x,node_y,sat_loc_x,sat_loc_y,phi,n,time_update,phi0)
A = ones(1,sat_num);
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
D = 550;
%satellite coordinates 
s = zeros(sat_num,length(t));
for n = 1:sat_num
    s(n,:) = A(n)*exp(1i*(w*n*time_update+phi(n)));
end

RNM = zeros(1,sat_num);
R = zeros(1,sat_num);
phi_moving = zeros(1,sat_num);
delta_phi = zeros(1,sat_num);
signal_receive = zeros(1,length(t));
%%  
        for m = 1:node_num
            for k = 1:sat_num
                R(k) = sqrt((y_ground(node_y(m))-sat_loc_y(k)).^2+(x_ground(node_x(m))-sat_loc_x(k)).^2+D^2);
                delta_phi(k) = R(k)*2*pi/lambda;
                signal_receive = signal_receive + s(k,:)*exp(1i*delta_phi(k));
            end
        cor(m) = signal_receive*signal_receive'/length(signal_receive);
        signal_receive = zeros(1,length(t));
        end
   power = mean(cor);
end

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

function satelliteplot(sattochoose_x,sattochoose_y,sat_x_record,sat_y_record,n)
    scatter(sattochoose_x,sattochoose_y,'x')
    hold on
    scatter(sat_x_record(n,:), sat_y_record(n,:),50,'b');
    legend('所有可选卫星','通过算法选取')
    pause(0.5)
    drawnow
    clf
end