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
node_x = x_ground(randi(length(x_ground)));
node_y = y_ground(randi(length(y_ground)));
H_sat = 10; %卫星分布范围
step_sat = H_sat/201;
sat_available = 10;
v_sat = 7.8;%卫星运动速度为7.8km/s
sat_x = -H_sat/2:step_sat:H_sat/2;
sat_y = -H_sat/2:step_sat:H_sat/2;
D = 550;%orbit 
nvars = 2*sat_num;%优化目标为相位和四颗卫星的标号
power_record = [];%用来记录每次更新的地面功率
outtime_record = [];%用来记录每次选星后的持续时间

figure%卫星图

for a = 1:5
% 产生可选卫星的坐标
sattochoose_x = sat_x(randperm(length(sat_x), sat_available));%可选卫星的坐标
sattochoose_y = sat_y(randperm(length(sat_y), sat_available));

% 产生可选卫星的轨道倾角（初始位置为随机生成的位置），用来计算卫星运动轨迹
sattrack_theta =  0 + rand(1,sat_available) * (90 - 0);%0到90度内的随机数
tic;
estimate_calculation_time = 0;%初始化计算时间/(以估计的运算时间及逆行计算)
% 使用 ga 函数进行单目标优化
targetFunction = @(x)singlenode_power(x,sat_num,node_x,node_y,sattochoose_x,sattochoose_y);
%为了避免重复，将选星分为四个区间，每个变量从一个区间中选取
UB =[pi*ones(1,sat_num)   round(sat_available/4) round(sat_available/2)      round(sat_available*3/4)    sat_available];
LB =[-pi*ones(1,sat_num)  1                      round(sat_available/4)+0.01 round(sat_available/2)+0.01 round(sat_available*3/4)+0.01];
options = optimoptions('ga','PopulationSize',50,'MaxGenerations',10);%default: ps 200; mg 200*number of variable
[x, fval,exitflag,output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,[],options);
calculation_time = toc

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
disp(['节点处功率值(NM)：',num2str(-fval)]) 
% [power_moving,sat_xmoving,sat_ymoving,phi_moving] = singlenode_power(x,sat_num,v_sat,calculation_time,sattrack_theta,node_x,node_y,x_sat,y_sat);
% disp(['节点处功率值(M)：',num2str(-power_moving)])
% disp(['卫星x方向位移',num2str(sat_xmoving)])
% disp(['卫星y方向位移',num2str(sat_ymoving)])
% disp(['卫星位移',num2str(calculation_time*v_sat)])
% disp(['相位变化',num2str(phi_moving)])

%%%%得出选择的四颗卫星和各自的相位后，开始根据卫星运动实时更新卫星权重
R0 = sqrt((node_x-optimalXloc).^2 + (node_y-optimalYloc).^2 + D^2);
sat_xmoving = v_sat*calculation_time*cosd(sattrack_theta(round(optimalidx)));
sat_ymoving = v_sat*calculation_time*sind(sattrack_theta(round(optimalidx)));
sat_xstart = optimalXloc + sat_xmoving;
sat_ystart = optimalYloc + sat_ymoving;%后续更新的初始卫星位置
R1 = sqrt((node_x-sat_xstart).^2 + (node_y-sat_ystart).^2 + D^2);%初始距离
phi_moving = 2*pi*(R1 - R0)/lambda;
phi = optimalPhi - phi_moving;%后续更新的初始相位

time_update = 0.01;%更新间隔
times = 100;%更新次数
sat_x_record = zeros(times,sat_num);
sat_y_record = zeros(times,sat_num);%用来记录四颗卫星的轨迹
phi_record = zeros(times,sat_num);%用来记录四颗卫星的相位

sat_x_record(1,:) = sat_xstart;
sat_y_record(1,:) = sat_ystart;
phi_record(1,:) = phi;
power_record = [power_record fval];

for n = 2:times
    sat_x_record(n,:) = sat_xstart + (n-1)*time_update*v_sat*cosd(sattrack_theta(round(optimalidx))); 
    sat_y_record(n,:) = sat_ystart + (n-1)*time_update*v_sat*sind(sattrack_theta(round(optimalidx))); 
    R_now = sqrt((node_x-sat_x_record(n,:)).^2 + (node_y-sat_y_record(n,:)).^2 + D^2);
    phi_record(n,:) = phi_record(n-1,:)-(R_now - R1)*2*pi/lambda;%根据之前的相位更新当前相位
    R1 = R_now;
    power_record = [power_record powerrecording(sat_num,node_x,node_y,sat_x_record(n,:),sat_y_record(n,:),phi_record(n,:))];
    x_out = any(sat_x_record(n,:) > H_sat);
    y_out = any(sat_y_record(n,:) > H_sat);
    %satelliteplot(sattochoose_x,sattochoose_y,sat_x_record,sat_y_record,n);
    
    if (x_out) || (y_out)
        break
    end
end
outtime_record = [outtime_record n];
end

powerplot_segmented(outtime_record,power_record);%绘制功率的分段曲线























%%%单节点功率计算
function [power,sat_xmoving,sat_ymoving,phi_moving] = singlenode_power(xx,sat_num,x,y,x_sat,y_sat)
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
t = 0:1e-10:1e-9;
s = zeros(sat_num,length(t));
for n = 1:sat_num
    s(n,:) = A(n)*exp(1i*(2*pi*frequency*t+phi(n)));
end

RNM = zeros(1,sat_num);
R = zeros(1,sat_num);
phi_moving = zeros(1,sat_num);
delta_phi = zeros(1,sat_num);
signal_receive = zeros(1,length(t));
%%
        for k = 1:sat_num 
            R(k) = sqrt((y-sat_loc_y(k)).^2+(x-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            %freespaceloss = lambda/(R(k)*4*pi);
            %归一化自由空间损耗,以轨道高度为基准
            freespaceloss = D/R(k);
            signal_receive = signal_receive + freespaceloss*s(k,:)*exp(1i*delta_phi(k));%考虑弗里斯传输方程给出的幅度衰减
        end
        power = -signal_receive*signal_receive'/length(signal_receive);%算法为最小化，故取负
        signal_receive = zeros(1,length(t));
end


%更新运动过程的power
function [power] = powerrecording(sat_num,x,y,sat_loc_x,sat_loc_y,phi)
A = ones(1,sat_num);
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
D = 550;
%satellite coordinates 
t = 0:1e-10:1e-9;
s = zeros(sat_num,length(t));
for n = 1:sat_num
    s(n,:) = A(n)*exp(1i*(2*pi*frequency*t+phi(n)));
end
RNM = zeros(1,sat_num);
R = zeros(1,sat_num);
phi_moving = zeros(1,sat_num);
delta_phi = zeros(1,sat_num);
signal_receive = zeros(1,length(t));
%%
        for k = 1:sat_num 
            R(k) = sqrt((y-sat_loc_y(k)).^2+(x-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            %freespaceloss = lambda/(R(k)*4*pi);
            %归一化自由空间损耗,以轨道高度为基准
            freespaceloss = D/R(k);
            signal_receive = signal_receive + freespaceloss*s(k,:)*exp(1i*delta_phi(k));%考虑弗里斯传输方程给出的幅度衰减
        end
        power = -signal_receive*signal_receive'/length(signal_receive);%算法为最小化，故取负
        signal_receive = zeros(1,length(t));
end

function powerplot_segmented(outtime_record,power_record)

% Initialize figure
figure;

startIndex = 1;  % Start index for each segment

for i = 1:length(outtime_record)
    endIndex = startIndex + outtime_record(i) - 1;  % End index for each segment
    
    % Extract segment from power_record
    segment_power = power_record(startIndex:endIndex);
    
    % Plot the segment
    plot(startIndex:endIndex, segment_power, 'LineWidth', 2);
    
    hold on;  % Hold the plot for the next segment
    
    startIndex = endIndex + 1;  % Update the start index for the next segment
end

hold off;  % Release the hold after all segments are plotted

title('Power Record Segments');
xlabel('Data Point');
ylabel('Power');
legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5');
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